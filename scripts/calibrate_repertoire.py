"""Bounded generation-only reconstruction; all trials are fresh native GenAIRR draws.

Run under profile_build.py. No reads, MS, Pgen, scoring, resume or feature caches.
"""
import argparse,csv,gzip,hashlib,json,math,os,random,shutil,subprocess,sys,time
from collections import Counter
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path


def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def inherited_cpu_budget(requested):
    """Use the resource monitor's validated allocation without scoring dependencies."""
    inherited=os.environ.get('POLARS_MAX_THREADS','')
    if not inherited.isascii() or not inherited.isdigit() or int(inherited)<1:
        raise ValueError('run under profile_build.py: positive inherited POLARS_MAX_THREADS budget required')
    if type(requested) is not int or requested<1:raise ValueError('cpus must be a positive integer')
    return min(requested,int(inherited))


def largest_remainder(counts,n):
    if type(n) is not int or n<0 or not counts or any(type(x) is not int or x<0 for x in counts.values()) or sum(counts.values())<=0:raise ValueError('invalid integer allocation')
    total=sum(counts.values());out={k:n*v//total for k,v in counts.items()}
    order=sorted(counts,key=lambda k:(-(n*counts[k]%total),k))
    for k in order[:n-sum(out.values())]:out[k]+=1
    return out


def wasserstein(a,b):
    if not a or not b or any(type(k) is not int or type(v) is not int or v<0 for h in [a,b] for k,v in h.items()) or min(sum(a.values()),sum(b.values()))<=0:raise ValueError('invalid integer histogram')
    na,nb=sum(a.values()),sum(b.values());ca=cb=num=0
    for x in range(min(min(a),min(b)),max(max(a),max(b))):
        ca+=a.get(x,0);cb+=b.get(x,0);num+=abs(ca*nb-cb*na)
    return dict(numerator=num,denominator=na*nb,distance_AA=num/(na*nb))


def validate_lambda_starts(starts):
    if (len(starts)!=3 or any(type(x) not in (int,float) or not 0<=x<=sys.float_info.max or not math.isfinite(x) for x in starts)
            or not starts[0]<starts[1]<starts[2] or starts[2]<=0):
        raise ValueError('lambda-starts require three finite strictly increasing nonnegative values')
    return tuple(float(x) for x in starts)


def next_lambda(results,starts=(0.,.1,.4)):
    starts=validate_lambda_starts(starts)
    if len(results)<3:return starts[len(results)]
    points=sorted(results,key=lambda r:r['lambda']);best=min(results,key=lambda r:(r['W1']['distance_AA'],r['lambda']))
    i=points.index(best);intervals=[]
    if i:intervals.append((points[i-1]['lambda'],best['lambda']))
    if i+1<len(points):intervals.append((best['lambda'],points[i+1]['lambda']))
    available=[];seen={Decimal(str(r['lambda'])) for r in results}
    for left,right in intervals:
        lo,hi=Decimal(str(left)),Decimal(str(right));midpoint=(lo+hi)/2
        rounded=midpoint.quantize(Decimal(1).scaleb(midpoint.adjusted()-1),rounding=ROUND_HALF_UP)
        if lo<rounded<hi and rounded not in seen:available.append((lo,hi,rounded))
    if not available:return None
    return float(min(available,key=lambda p:(-(p[1]-p[0]),p[0]))[2])


def tilted_parameters(initial,lam):
    """Keep the authored positive measure; order small masses first for native CDF safety."""
    if type(lam) not in (int,float) or not math.isfinite(lam) or lam<0:raise ValueError('invalid lambda')
    params=dict(initial)
    for key in ['np1_lengths','np2_lengths']:
        weights=[w*math.exp(-lam*n) for n,w in initial[key]];total=sum(weights)
        pairs=[[n,w/total] for (n,_),w in zip(initial[key],weights)]
        if any(not math.isfinite(w) or w<=0 for n,w in pairs):raise ValueError('NP tilt lost positive support')
        params[key]=sorted(pairs,key=lambda p:(p[1],p[0]))
    return params


def preflight_plans(initial,starts):
    """Compile every declared start before output creation; never sample sequences."""
    from GenAIRR import Experiment
    receipts=[]
    for lam in validate_lambda_starts(starts):
        params=tilted_parameters(initial,lam)
        plan=Experiment.on('human_igh').recombine(**params).productive_only().compile()
        receipts.append(dict(**{'lambda':lam},status='compiled',pass_plan=list(plan.pass_plan),
                             parameter_sha256=hashlib.sha256(json.dumps(params,sort_keys=True,separators=(',',':')).encode()).hexdigest()))
    return receipts


def targets(path,condition,population):
    with Path(path).open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
    selected=[r for r in rows if r['condition']==condition and r['population']==population]
    out={}
    for q in ['V_gene_call','J_gene_call','C_gene_call','anchor_inclusive_junction_AA_length','supplied_annotation_variable_AA_length']:
        hist={}
        for r in selected:
            if r['quantity']!=q:continue
            key=r['category'];n=int(r['numerator'])
            if key in hist or n<0 or int(r['denominator'])!=33040 or int(r['n_observations'])!=33040:raise ValueError('invalid target bin/denominator')
            hist[key]=n
        if sum(hist.values())!=33040:raise ValueError('target histogram must retain all33040 records: '+q)
        out[q]=hist
    length_bins=out['anchor_inclusive_junction_AA_length']
    if length_bins.get('missing',0):raise ValueError('target junction lengths missing')
    lengths={int(k):v for k,v in length_bins.items() if k!='missing'}
    if any(k<1 for k in lengths):raise ValueError('invalid target junction length')
    if set(out['C_gene_call'])!={'IGHM','IGHG1','IGHG3','IGHA1','IGHE','missing'}:raise ValueError('unsupported target constant bins')
    return out,lengths


def tv(observed,target):
    keys=sorted(set(observed)|set(target));no,nt=sum(observed.values()),sum(target.values())
    rows=[dict(category=k,observed=observed.get(k,0),target=target.get(k,0),observed_denominator=no,target_denominator=nt) for k in keys]
    return dict(value=sum(abs(observed.get(k,0)/no-target.get(k,0)/nt) for k in keys)/2 if no and nt else None,bins=rows)


def constant_inputs(target,seed,layout_ids=None):
    if layout_ids is None:layout_ids=[f'family{f}_node{n}' for f in range(800) for n in range(1 if f<600 else 3)]
    if (not isinstance(layout_ids,(list,tuple)) or not layout_ids or any(not isinstance(x,str) or not x for x in layout_ids) or len(set(layout_ids))!=len(layout_ids)):
        raise ValueError('layout IDs must be nonempty unique antibody strings')
    labels=largest_remainder(target,len(layout_ids));known={k:v for k,v in target.items() if k!='missing'}
    imputed=largest_remainder(known,labels.get('missing',0));values=[]
    for gene,n in sorted(labels.items()):
        if gene!='missing':values.extend((gene,gene,False) for _ in range(n))
    for gene,n in sorted(imputed.items()):values.extend(('missing',gene,True) for _ in range(n))
    random.Random(seed).shuffle(values)
    return {name:v[1] for name,v in zip(layout_ids,values)},[dict(antibody_id=name,observed_C=v[0],physical_C=v[1],imputed=v[2]) for name,v in zip(layout_ids,values)]


def native_counter_summary(attempts):
    """Sum actual nonnegative integer event counts; report absent/invalid fields."""
    counters=[a.get('native_counters',{}) for a in attempts]
    if any(not isinstance(c,dict) for c in counters):raise ValueError('native_counters must be objects')
    result={}
    for key in sorted({k for c in counters for k in c if k.startswith('n_')}):
        valid=[c[key] for c in counters if key in c and type(c[key]) is int and c[key]>=0]
        missing=sum(key not in c or c[key] is None for c in counters)
        result[key]=dict(sum=sum(valid) if valid else None,n_valid=len(valid),n_missing=missing,n_invalid=len(counters)-len(valid)-missing)
    return result


def write_json(path,value):path.write_text(json.dumps(value,indent=2)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['targets','initial-model','reference-dir','additional-reference-dir','output']:p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--condition',required=True);p.add_argument('--population',required=True,choices=['background'])
    p.add_argument('--cpus',type=int,default=16)
    p.add_argument('--lambda-starts',type=float,nargs=3,default=[0.,.1,.4]);a=p.parse_args()
    starts=validate_lambda_starts(a.lambda_starts)
    repo=Path(__file__).resolve().parents[1]
    budget=inherited_cpu_budget(a.cpus)
    from generate_population import load_recombination_parameters, CONSTANT_MANIFEST
    from smoke_generate import reference_hashes
    initial,_=load_recombination_parameters(a.initial_model);model=json.loads(a.initial_model.read_text())
    target,length_target=targets(a.targets,a.condition,a.population)
    if a.output.exists():p.error('fresh absent output required')
    initial_plan_preflight=preflight_plans(initial,starts)
    inputs=[a.targets,a.initial_model,repo/'envs/generation/identity.json',*sorted((repo/'envs/generation').glob('*.json')),
            *[Path(f) for f in reference_hashes(a.reference_dir,all_loci=True)],
            *[a.additional_reference_dir/name for name in json.loads(CONSTANT_MANIFEST.read_text())['downloaded_records']],
            *[repo/'scripts'/n for n in ['calibrate_repertoire.py','generate_population.py','smoke_generate.py','validate_constant_candidates.py','genairr_gene_map.py']]]
    from importlib.resources import files
    package=Path(str(files('GenAIRR')))
    inputs.extend(f for f in package.rglob('*') if f.is_file() and '__pycache__' not in f.parts and f.suffix!='.pyc')
    inputs.extend((repo/'envs/generation'/name) for name in ['requirements.txt','requirements.in'])
    pins={str(f.resolve()):sha(f) for f in inputs};a.output.mkdir(parents=True)
    write_json(a.output/'target_histograms.json',target)
    allocation,labels=constant_inputs(target['C_gene_call'],20261111)
    write_json(a.output/'allocation.json',dict(schema_version=1,allocation=allocation))
    with (a.output/'C_observation_labels.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(labels[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(labels)
    results=[];start=time.monotonic();receipt=dict(status='running',np_ordering_policy='ascending(weight,length); no floor or additional rescaling; seeded inverse-CDF realizations change with ordering',initial_plan_preflight=initial_plan_preflight,lambda_starts=list(starts),lambda_bounds=[starts[0],starts[-1]],max_trials=8,stop_reason=None,source_input_sha256=pins,cpu_budget=budget,seed=20261111,retention='gzip full final truth; parameters, labels and summaries; remove only trial-generated FASTAs and copied references after hashing',trials=results)
    def save():write_json(a.output/'execution.json',receipt)
    save()
    try:
        stop_reason='max_trials'
        for i in range(8):
            if time.monotonic()-start>=14400:raise RuntimeError('aggregate four-hour bound')
            lam=next_lambda(results,starts)
            if lam is None:
                stop_reason='parameter_precision_limit';break
            trial=a.output/f'trial{i+1:02d}';trial.mkdir();modified=json.loads(json.dumps(model))
            modified['parameters']=tilted_parameters(initial,lam)
            write_json(trial/'model.json',modified);dest=trial/'native';tick=time.monotonic()
            command=[sys.executable,str(repo/'scripts/generate_population.py'),'--output',str(dest),'--reference-dir',str(a.reference_dir),'--additional-reference-dir',str(a.additional_reference_dir),'--total-antibodies','1200','--memory-records','600','--seed','20261111','--mutation-rate','0.01','--pairs-per-chain','100','--recombination-parameters',str(trial/'model.json'),'--heavy-constant-allocation',str(a.output/'allocation.json')]
            with (trial/'generator.log').open('w') as log:subprocess.run(command,check=True,stdout=log,stderr=subprocess.STDOUT,timeout=min(1800,14400-(time.monotonic()-start)))
            generation_wall_seconds=time.monotonic()-tick
            truth_path=dest/'truth.json';truth=json.loads(truth_path.read_text());rows=truth['chains'];heavy=[r for r in rows if r['locus']=='IGH']
            if len(rows)!=2400 or len(heavy)!=1200 or len({r['chain_id'] for r in rows})!=2400 or len(truth['latent_ancestors'])!=400:raise ValueError('native population invariant')
            if Counter(r['locus'] for r in rows)!=Counter(IGH=1200,IGK=600,IGL=600):raise ValueError('locus counts differ')
            pairs={}
            for r in rows:pairs.setdefault(r['antibody_id'],[]).append(r['locus'])
            if len(pairs)!=1200 or any(len(v)!=2 or 'IGH' not in v or len(set(v))!=2 for v in pairs.values()):raise ValueError('pairing invariant')
            lengths=Counter();vg=Counter();jg=Counter();vl=Counter();vnt=Counter();jnt=Counter()
            for r in heavy:
                junction=r['airr_truth']['junction_aa']
                if not junction or junction[0]!='C' or junction[-1] not in 'FW':raise ValueError('native anchored junction missing')
                vnt[r['airr_truth']['v_sequence_end']-r['airr_truth']['v_sequence_start']]+=1;jnt[r['airr_truth']['j_sequence_end']-r['airr_truth']['j_sequence_start']]+=1
                lengths[len(junction)]+=1;vg[r['canonical_v_gene'] or 'unresolved:'+r['canonical_v_match_status']]+=1;jg[r['airr_truth']['j_call'].split('*')[0]]+=1;vl[len(r['variable_aa'])]+=1
                if r['constant_nt_gene']!=allocation[r['antibody_id']]:raise ValueError('constant allocation drift')
            hist=dict(V_nt_length=dict(vnt),J_nt_length=dict(jnt),junction_length=dict(lengths),variable_length=dict(vl),V=dict(vg),J=dict(jg),locus=dict(Counter(r['locus'] for r in rows)),C_physical=dict(Counter(r['constant_nt_gene'] for r in heavy)),C_projected=dict(Counter(r['observed_C'] for r in labels)),C_imputed=dict(Counter(r['physical_C'] for r in labels if r['imputed'])),components=dict(Counter(r['component'] for r in heavy)),V_mapping_status=dict(Counter(r['canonical_v_match_status'] for r in heavy)),family_size=dict(Counter(Counter(r['family_id'] for r in heavy).values())))
            for field in ['C_physical','C_imputed']:
                for gene in target['C_gene_call']:
                    if gene!='missing':hist[field].setdefault(gene,0)
            for gene in target['C_gene_call']:hist['C_projected'].setdefault(gene,0)
            metrics={}
            for key,generated,q in [('V',vg,'V_gene_call'),('J',jg,'J_gene_call')]:
                observed=target[q];supported=set(observed)-{'missing'}
                metrics[key]=dict(full=tv(generated,observed),conditional_observed_support=tv({k:v for k,v in generated.items() if k in supported},{k:v for k,v in observed.items() if k in supported}),generated_residual_records=sum(v for k,v in generated.items() if k not in supported),target_records_in_zero_generated_bins=sum(v for k,v in observed.items() if k not in generated))
            item=dict(trial=i+1,**{'lambda':lam},W1=wasserstein(dict(lengths),length_target),generation_wall_seconds=generation_wall_seconds,command=command,latent_chains=len(truth['latent_ancestors']),histograms=hist,gene_diagnostics=metrics,rejected_draws=truth['rejected_draws'],rejection_reasons=dict(Counter(str(x['rejection_reason']) for x in truth['attempts'])),SHM={l:dict(chains=len(rr),substitutions=sum(r['root_distance']['substitutions'] for r in rr),denominator_nt=sum(r['root_distance']['denominator_nt'] for r in rr)) for l in ['IGH','IGK','IGL'] for rr in [[r for r in rows if r['locus']==l]]})
            item['native_integer_counts_over_recorded_attempts']=native_counter_summary(truth['attempts'])
            if Counter(r['component'] for r in heavy)!=Counter(unmutated_singleton=600,related_shm_family=600):raise ValueError('component count differs')
            raw_hash=sha(truth_path);raw_size=truth_path.stat().st_size;packed=dest/'truth.json.gz'
            with truth_path.open('rb') as source,packed.open('wb') as output:
                with gzip.GzipFile(filename='',fileobj=output,mode='wb',mtime=0) as sink:shutil.copyfileobj(source,sink)
            with gzip.open(packed,'rb') as source:
                if hashlib.file_digest(source,'sha256').hexdigest()!=raw_hash:raise ValueError('gzip roundtrip mismatch')
            item['truth']=dict(raw_sha256=raw_hash,raw_bytes=raw_size,gzip_sha256=sha(packed),gzip_bytes=packed.stat().st_size);truth_path.unlink()
            removed={str(f.relative_to(trial)):sha(f) for f in [*dest.glob('*.fasta'),*[f for f in (dest/'references').rglob('*') if f.is_file()]]}
            for relative,digest in removed.items():
                f=trial/relative
                if sha(f)!=digest:raise ValueError('generated artifact changed before cleanup')
                if f.is_relative_to(dest/'references'):
                    key=str(f.relative_to(dest/'references'))
                    if key in truth['reference_sha256'] and digest!=truth['reference_sha256'][key]:raise ValueError('copied reference hash mismatch')
            for f in dest.glob('*.fasta'):f.unlink()
            shutil.rmtree(dest/'references');item['removed_generated_artifact_sha256']=removed
            item['wall_seconds']=time.monotonic()-tick
            write_json(trial/'summary.json',item);results.append(item);save();print(json.dumps(dict(trial=i+1,**{'lambda':lam},W1=item['W1'])),flush=True)
        if any(sha(Path(f))!=h for f,h in pins.items()):raise RuntimeError('source/input drift')
        best=min(results,key=lambda r:(r['W1']['distance_AA'],r['lambda'],r['trial']));receipt.update(status='complete',stop_reason=stop_reason,selected_trial=best['trial'],selected_lambda=best['lambda'],selected_at_boundary=best['lambda'] in (starts[0],starts[-1]),wall_seconds=time.monotonic()-start)
    except BaseException as error:receipt.update(status='failed',error=repr(error));raise
    finally:
        receipt['artifact_sha256']={str(f.relative_to(a.output)):sha(f) for f in a.output.rglob('*') if f.is_file() and f.name!='execution.json'};save()


if __name__=='__main__':main()
