"""Require explicit background roles before the legacy supplied-panel adapter."""
from simulate_airr_ms import main, parse_args

if __name__ == '__main__':
    args = parse_args()
    import vdjaqp
    vdjaqp.resolve_cpu_budget(args.cpus)
    import polars as pl
    candidates = pl.read_parquet(args.repertoire)
    if ('background_member' not in candidates.columns
            or candidates.schema['background_member'] != pl.Boolean
            or candidates['background_member'].null_count()):
        raise ValueError('workflow requires explicit Boolean background_member')
    panel = pl.read_csv(args.target_ids)['clone_id']
    if candidates.filter(pl.col('clone_id').is_in(panel.implode()))['background_member'].any():
        raise ValueError('panel must contain premarked computational references, never observed background targets')
    main()
    from revalidate_airr_background_membership import scoped_metrics
    positives=set(pl.read_parquet(args.output/'physical_truth.parquet')['chain_id'])
    references=set(pl.read_parquet(args.output/'pooled_panel_coverage.parquet')['chain_id'])
    if len(references)!=len(panel):raise ValueError('canonical exported panel membership changed')
    pl.DataFrame([dict(enzyme=enzyme,**row)
                  for enzyme in ['Trypsin','Chymotrypsin','pooled']
                  for row in scoped_metrics(pl.read_parquet(args.output/(enzyme+'_scores.parquet')),positives,references)]).write_csv(args.output/'classification_statistics.tsv',separator='\t')
    from compress_airr_ms_artifacts import compress_case
    compress_case(args.output)
