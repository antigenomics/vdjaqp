"""Export native sequence-view evidence for the declared full candidate population.

Retained-context mode formats this fresh DAG's placements and scores without
rebuilding evidence. The historical default builds Evidence and selects references
and four high-scoring background records; it is not the primary raw-science route.
Display selection is neither population sampling nor a new acceptance threshold.
"""
import argparse
import hashlib
import json
from pathlib import Path


def export_retained(input_root, context_root, display_table, output,
                    expected_candidates, case='gel', cpus=None):
    """Format this run's retained placements/scores; never rebuild Evidence."""
    import vdjaqp
    budget=vdjaqp.resolve_cpu_budget(cpus)
    import numpy as np
    import polars as pl
    from vdjaqp.report import placement_depth
    import inspect
    if output.exists():
        raise FileExistsError('fresh display output required')
    def sha(path):
        with path.open('rb') as stream:
            return hashlib.file_digest(stream,'sha256').hexdigest()
    analysis=context_root/'analysis'
    provenance_path=analysis/'provenance.json'
    provenance_sha=sha(provenance_path)
    provenance=json.loads(provenance_path.read_text())
    if provenance['status']!='complete':
        raise ValueError('complete context evidence is required')
    entries=[row for row in provenance['cases'] if row['case']==case]
    if len(entries)!=1:
        raise ValueError('one declared context case is required')
    mapping=context_root/'peptide_mapping'
    receipt=mapping/'mapping_receipt.json'
    receipt_sha=sha(receipt)
    if receipt_sha!=entries[0]['mapping_receipt_sha256']:
        raise ValueError('mapping receipt differs from scored context')
    mapped=json.loads(receipt.read_text())
    if mapped['status']!='complete' or mapped['query_role']!='identity_eligible':
        raise ValueError('identity-query placement receipt required')
    sources={str(path):sha(path) for path in [Path(__file__).resolve(),Path(inspect.getsourcefile(placement_depth)).resolve()]}
    sources.update({str(provenance_path):provenance_sha,str(receipt):receipt_sha})
    for name,expected in mapped['artifact_sha256'].items():
        path=mapping/name
        if sha(path)!=expected:
            raise ValueError(f'placement artifact drift: {name}')
        sources[str(path)]=expected
    candidate_path=analysis/case/'candidates.parquet'
    score_path=analysis/case/'alignment_union_scores.parquet'
    for path in (candidate_path,score_path):
        expected=provenance['artifact_sha256'][str(path.relative_to(analysis))]
        if sha(path)!=expected:
            raise ValueError(f'analysis artifact drift: {path.name}')
        sources[str(path)]=expected
    targets=input_root/'ground_truth/target_clonotypes.csv'
    target_sha=sha(targets)
    if target_sha!=provenance['raw_sha256']['ground_truth/target_clonotypes.csv']:
        raise ValueError('reference metadata drift')
    sources[str(targets)]=target_sha
    refs=set(pl.read_csv(targets)['clone_id'])
    sources[str(display_table)]=sha(display_table)
    selected=pl.read_csv(display_table,separator='\t',schema_overrides={'clone_id':pl.String,'display_id':pl.String})
    if selected.height==0 or selected.select('clone_id','display_id').null_count().sum_horizontal().item() or any(selected[c].n_unique()!=selected.height for c in ['clone_id','display_id']):
        raise ValueError('nonempty unique display/query IDs are required')
    clones=pl.read_parquet(candidate_path)
    if clones.height!=expected_candidates or clones['clone_id'].n_unique()!=expected_candidates:
        raise ValueError('declared full candidate population differs')
    queries=pl.read_parquet(mapping/'queries.parquet')
    if not clones.select('clone_id','seq').sort('clone_id').equals(queries.sort('clone_id')):
        raise ValueError('mapping/query sequence mismatch')
    chosen=selected.join(clones,on='clone_id',how='left',validate='1:1',maintain_order='left')
    if chosen['seq'].null_count():
        raise ValueError('display table names an unavailable query')
    scores=pl.read_parquet(score_path,columns=['clone_id','confirmed_nats','junction_nats','framework_only'])
    if scores.height!=expected_candidates or scores['clone_id'].n_unique()!=expected_candidates:
        raise ValueError('complete candidate scores required')
    hits=pl.read_parquet(mapping/'placements.parquet')
    ambiguity=hits.group_by('peptide_id').agg(pl.col('clone_id').n_unique().alias('compatible_candidate_records'))
    observations=pl.read_parquet(mapping/'observations.parquet')
    hits=hits.filter(pl.col('clone_id').is_in(selected['clone_id'].implode())).join(
        observations.rename({'seq':'peptide'}),on='peptide_id',how='left',validate='m:1').join(ambiguity,on='peptide_id',validate='m:1')
    if hits['peptide'].null_count():
        raise ValueError('retained placement lacks its observation sequence')
    table=[];residues=[];placements=[]
    for row in chosen.iter_rows(named=True):
        label=row['display_id'];seq=row['seq'];start=int(row['cdr3_start']);end=int(row['cdr3_end'])
        if not 0<=start<end<=len(seq):
            raise ValueError('display query requires declared complete junction geometry')
        hh=hits.filter(pl.col('clone_id')==row['clone_id'])
        if hh.filter((pl.col('start')<0)|(pl.col('end')>len(seq))|((pl.col('end')-pl.col('start'))!=pl.col('peptide').str.len_chars())).height:
            raise ValueError('invalid ungapped native placement geometry')
        depth,agreeing=placement_depth(seq,hh)
        covered=agreeing>0
        support=row.get('nt_support')
        if support is None:
            support=np.asarray(row['region_ids'])!=8
            support_kind='supplied-region proxy, not nucleotide ancestry'
        else:
            support=np.asarray(support)
            if support.dtype!=np.bool_ or support.shape!=(len(seq),):
                raise ValueError('one Boolean nucleotide-support value per residue required')
            support_kind='explicit supplied per-residue support'
        table.append(dict(display_id=label,clone_id=row['clone_id'],is_reference=row['clone_id'] in refs,
            defined_aa=len(seq),covered_aa=int(covered.sum()),junction_defined_aa=end-start,
            junction_covered_aa=int(covered[start:end].sum()),compatible_placements=hh.height,
            compatible_peptide_rows=hh['peptide_id'].n_unique(),sequence_origin=row.get('sequence_origin'),
            coverage_kind=row.get('coverage_kind'),junction_core_status=row.get('junction_core_status'),
            nt_support_provenance=support_kind))
        for at,aa in enumerate(seq):
            residues.append(dict(display_id=label,position=at,junction_relative_position=at-start,
                aa=aa,region_id=int(row['region_ids'][at]),covered=bool(covered[at]),
                nt_supported=bool(support[at]),agreeing_depth=int(agreeing[at]),
                compatible_depth=int(depth[at]),in_junction=start<=at<end))
        placements.append(hh.with_columns(pl.lit(label).alias('display_id')))
    table=pl.DataFrame(table).join(scores,on='clone_id',validate='1:1',maintain_order='left')
    if table['confirmed_nats'].null_count():
        raise ValueError('display query has no retained score')
    if any(sha(Path(path))!=expected for path,expected in sources.items()):
        raise ValueError('retained input drift while formatting')
    output.mkdir(parents=True)
    table.write_csv(output/'clonotypes.tsv',separator='\t')
    pl.DataFrame(residues).write_csv(output/'residues.tsv',separator='\t')
    pl.concat(placements,how='diagonal_relaxed').write_parquet(output/'placements.parquet')
    observations.write_parquet(output/'peptides.parquet')
    (output/'manifest.json').write_text(json.dumps(dict(status='complete',case=case,cpus=budget,
        search_or_scoring_repeated=False,candidate_count=expected_candidates,selected_count=selected.height,
        post_analysis='coverage from retained placements using vdjaqp.report.placement_depth; retained-score join; declared display subset',
        selection='explicit declared display table; not a population sample or accepted-identity set',
        coverage='agreeing compatible placements under I/L equivalence; no unique credit implied',
        coordinates='zero-based half-open; junction includes conserved anchors; residue positions are junction-start relative',
        ambiguity_universe='complete identity-query candidate records; partial compatibility-only intervals excluded',
        source_sha256=sources,artifact_sha256={path.name:sha(path) for path in output.iterdir() if path.is_file()}),indent=2)+'\n')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--expected-candidates', required=True, type=int)
    p.add_argument('--cpus', type=int)
    p.add_argument('--retained-context',type=Path,help='Format existing placements/scores from this fresh DAG, without a new Evidence build')
    p.add_argument('--display-table',type=Path,help='Declared TSV with unique clone_id/display_id for retained-context mode')
    p.add_argument('--case',default='gel')
    a = p.parse_args()
    if a.expected_candidates < 1:
        p.error('expected candidate count must be positive')
    if a.output.exists():
        p.error('fresh output required')
    if a.retained_context is not None:
        if a.display_table is None:
            p.error('retained context requires an explicit display table')
        export_retained(a.input_root,a.retained_context,a.display_table,a.output,a.expected_candidates,a.case,a.cpus)
        return
    import vdjaqp
    budget = vdjaqp.resolve_cpu_budget(a.cpus)
    import numpy as np
    import polars as pl
    from vdjaqp.data import load_targets
    from vdjaqp.pipeline import Evidence
    from vdjaqp.report import profile
    files = [q for q in a.input_root.rglob('*') if q.is_file()]
    sha = lambda q: hashlib.sha256(q.read_bytes()).hexdigest()
    raw = {str(q.relative_to(a.input_root)): sha(q) for q in files}
    ev = Evidence.build('dataset_mixa_gel', root=a.input_root)
    if ev.clonotypes.height != a.expected_candidates:
        raise ValueError('declared population differs')
    refs = set(load_targets(a.input_root)['clone_id'])
    scored = ev.score().sort(['score', 'clone_id'], descending=[True, False])
    selected = sorted(refs) + scored.filter(~pl.col('clone_id').is_in(refs))['clone_id'].head(4).to_list()
    residues, placements, table = [], [], []
    meta = {r['clone_id']: r for r in ev.clonotypes.filter(pl.col('clone_id').is_in(selected)).iter_rows(named=True)}
    for i, clone_id in enumerate(selected):
        label = f'Ref{i+1:02d}' if clone_id in refs else f'Bio{i-len(refs)+1:02d}'
        prof = profile(clone_id, ev, targets=refs)
        m = meta[clone_id]
        start, end = int(m['cdr3_start']), int(m['cdr3_end'])
        cov = prof.agreeing_depth > 0
        table.append(dict(display_id=label, clone_id=clone_id, is_reference=clone_id in refs,
                          defined_aa=len(prof.seq), covered_aa=int(cov.sum()),
                          junction_defined_aa=end-start, junction_covered_aa=int(cov[start:end].sum()),
                          compatible_placements=prof.peptides.height, compatible_peptide_rows=prof.peptides['peptide_id'].n_unique(),
                          sequence_origin=m.get('sequence_origin'),coverage_kind=m.get('coverage_kind'),junction_core_status=m.get('junction_core_status')))
        for j, aa in enumerate(prof.seq):
            residues.append(dict(display_id=label, position=j, junction_relative_position=j-start,
                                 aa=aa, region_id=int(prof.regions[j]), covered=bool(cov[j]),nt_supported=bool(prof.nt_support[j]),
                                 agreeing_depth=int(prof.agreeing_depth[j]),
                                 compatible_depth=int(prof.depth[j]), in_junction=start<=j<end))
        if prof.peptides.height:
            placements.append(prof.peptides.with_columns(pl.lit(label).alias('display_id')))
    a.output.mkdir(parents=True)
    pl.DataFrame(table).write_csv(a.output/'clonotypes.tsv', separator='\t')
    pl.DataFrame(residues).write_csv(a.output/'residues.tsv', separator='\t')
    if placements:
        pl.concat(placements, how='diagonal_relaxed').write_parquet(a.output/'placements.parquet')
    ev.peptides.write_parquet(a.output/'peptides.parquet')
    if any(sha(a.input_root/n) != v for n, v in raw.items()):
        raise ValueError('input drift')
    (a.output/'manifest.json').write_text(json.dumps(dict(status='complete', cpus=budget,
        candidate_count=ev.clonotypes.height, selected_count=len(selected), raw_sha256=raw,
        selection='all 12 inserted references plus four highest full-score biological candidates; deterministic ties',
        coordinates='zero-based half-open; junction includes anchors; residue profile alignment is junction-start relative',
        coverage='at least one agreeing native peptide placement under I/L equivalence; compatible, not uniquely credited',
        artifact_sha256={q.name:sha(q) for q in a.output.iterdir() if q.is_file()}), indent=2)+'\n')


if __name__ == '__main__':
    main()
