"""Prefix validated mature model transcripts; no genotype draws, reads or annotation.

Run in the pinned generation environment under the outer resource monitor.
Coordinates are truth-only provenance and never candidate annotation boundaries.
"""
import argparse,csv,gzip,hashlib,importlib.util,json,shutil,time
from collections import Counter
from pathlib import Path
from Bio.Seq import translate

REPO=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('leader_context',REPO/'vdjaqp/transcript_context.py')
context_api=importlib.util.module_from_spec(spec);spec.loader.exec_module(context_api)


def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def load_contexts(path,expected_sha,source_path):
    if sha(path)!=expected_sha:raise ValueError('pinned context model hash mismatch')
    audit=json.loads((path.parent/'audit.json').read_text())
    if sha(source_path)!=audit['model_source_sha256']:raise ValueError('context source model mismatch')
    source=json.loads(source_path.read_text());native={r['native_allele']:r for r in source['native']};sources={r['source_id']:r for r in source['sources']}
    with path.open() as stream:rows=list(csv.DictReader(stream,delimiter='\t'))
    contexts={r['native_allele']:r for r in rows}
    if len(contexts)!=len(rows) or set(contexts)!=set(native) or len(contexts)!=audit['native_alleles']:raise ValueError('complete native context population required')
    if dict(Counter(r['method'] for r in rows))!=audit['final_methods']:raise ValueError('context method counts mismatch')
    for allele,r in contexts.items():
        donor=sources[r['selected_source_id']]
        if (r['method'] not in ['STRICT_EXACT','IMPUTED_GROUP','IMPUTED_ALIGNMENT'] or r['locus']!=native[allele]['locus'] or donor['locus']!=r['locus'] or donor['status']!='valid' or r['prefix_nt']!=donor['prefix_nt'] or r['leader_nt']!=donor['leader_nt'] or hashlib.sha256(r['prefix_nt'].encode()).hexdigest()!=r['prefix_sha256']):raise ValueError('context source/prefix mismatch')
    return contexts


def convert_chain(row,contexts):
    """Use only mature chain fields and native V truth needed for context lookup."""
    nt,mature,protein=row['variable_nt'],row['rnaseq_nt'],row['protein_aa']
    if not nt or not mature.startswith(nt) or len(mature)<=len(nt):raise ValueError('mature variable prefix mismatch')
    native=row['airr_truth'];allele=native['v_call'];context=contexts[allele]
    template=context_api.leader_rna_template(nt,mature[len(nt):],protein,context,v_trim_5=native['v_trim_5'],translate_nt=lambda s:str(translate(s)))
    provenance=dict(chain_id=row['chain_id'],native_v_allele=allele,locus=context['locus'],mature_nt_start=template['mature_nt_start'],leader_nt_start=template['leader_nt_start'],mature_nt_length=template['mature_nt_length'],method=context['method'],confidence=context['confidence'],source_id=context['selected_source_id'],prefix_sha256=context['prefix_sha256'])
    return template['rnaseq_nt'],provenance


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['truth-gzip','validation-receipt','contexts','context-source-model','output']:p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--contexts-sha256',required=True);a=p.parse_args();start=time.monotonic()
    if a.output.exists():p.error('fresh absent output required')
    validation=json.loads(a.validation_receipt.read_text());native=a.truth_gzip.parent
    if validation['status']!='complete' or validation['observed_chains']!=66080 or validation['background_antibodies']!=33040:raise ValueError('complete full-cohort validation required')
    if sha(a.truth_gzip)!=validation['truth']['gzip_sha256']:raise ValueError('native truth gzip hash mismatch')
    contexts=load_contexts(a.contexts,a.contexts_sha256,a.context_source_model)
    reference_pins=json.loads(a.context_source_model.read_text())['native_reference_sha256']
    for name,digest in reference_pins.items():
        matches=[v for path,v in validation['source_input_sha256'].items() if path.endswith('/builtin_dataconfigs/'+name)]
        if matches!=[digest]:raise ValueError('context model/native generator reference mismatch: '+name)
    paths=[a.truth_gzip,a.validation_receipt,a.contexts,a.context_source_model,a.contexts.parent/'audit.json',Path(__file__),REPO/'vdjaqp/transcript_context.py',*[native/name for name in ['read_counts.tsv','rnaseq_templates.fasta','amplicon_templates.fasta','proteins.fasta']]]
    pins={str(f.resolve()):sha(f) for f in paths}
    for name in ['read_counts.tsv','rnaseq_templates.fasta','amplicon_templates.fasta','proteins.fasta']:
        f=native/name
        if pins[str(f.resolve())]!=validation['source_input_sha256'][str(f.resolve())]:raise ValueError('validated native substrate changed: '+name)
    # Read native JSON once; raw digest is computed on the exact decompressed bytes.
    with gzip.open(a.truth_gzip,'rb') as stream:payload=stream.read()
    if hashlib.sha256(payload).hexdigest()!=validation['truth']['raw_sha256'] or len(payload)!=validation['truth']['raw_bytes']:raise ValueError('native truth roundtrip mismatch')
    truth=json.loads(payload);del payload
    rows=truth['chains'];ids=set();loci=Counter();methods=Counter()
    a.output.mkdir(parents=True)
    fields=['chain_id','native_v_allele','locus','mature_nt_start','leader_nt_start','mature_nt_length','method','confidence','source_id','prefix_sha256']
    with (a.output/'rnaseq_templates.fasta').open('w') as fasta,(a.output/'context_provenance.tsv').open('w') as table:
        writer=csv.DictWriter(table,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for row in rows:
            seq,provenance=convert_chain(row,contexts);cid=provenance['chain_id']
            if cid in ids:raise ValueError('duplicate native chain ID')
            ids.add(cid);loci[provenance['locus']]+=1;methods[provenance['method']]+=1
            fasta.write(f'>{cid}\n{seq}\n');writer.writerow(provenance)
    if len(ids)!=66080 or loci!=Counter(IGH=33040,IGK=16520,IGL=16520):raise ValueError('full native chain/locus population differs')
    seen=set()
    with (native/'read_counts.tsv').open() as stream:
        for cid,n in csv.reader(stream,delimiter='\t'):
            if cid in seen or cid not in ids or int(n)!=200:raise ValueError('prescribed native read allocation differs')
            seen.add(cid)
    if seen!=ids:raise ValueError('read model membership differs')
    shutil.copyfile(native/'read_counts.tsv',a.output/'read_counts.tsv')
    if sha(a.output/'read_counts.tsv')!=pins[str((native/'read_counts.tsv').resolve())]:raise ValueError('read model copy mismatch')
    if any(sha(Path(f))!=digest for f,digest in pins.items()):raise ValueError('source/model/native substrate changed')
    receipt=dict(status='complete',model='leader_imputed_v1',observed_chains=66080,locus_counts=dict(loci),context_methods=dict(methods),truth_only_coordinates=True,candidate_boundary_policy='numeric mature/leader offsets are truth provenance only; never feed candidate annotation boundaries',native_substrate_policy='all original mature chains, coding-only RNA/amp/proteins/truth bytes preserved; no redraw',read_counts_scope='copied prescribed100pairs per chain, not observed RNA abundance',source_input_sha256=pins,artifact_sha256={f.name:sha(f) for f in a.output.iterdir() if f.is_file()},wall_seconds=time.monotonic()-start)
    (a.output/'manifest.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:receipt[k] for k in ['status','observed_chains','locus_counts','context_methods','wall_seconds']}))


if __name__=='__main__':main()
