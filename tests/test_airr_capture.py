import numpy as np
import polars as pl
import pytest
from types import SimpleNamespace
from vdjaqp.airr_capture import capture_tables


def fixture():
    c=pl.DataFrame(dict(clone_id=['a','b','ref'],seq=['ACAWFG','ACAWFG','ACAWFG'],cdr3_seq=['CAWF']*3,cdr3_start=[1]*3,cdr3_end=[5]*3,region_ids=[[0,5,6,6,7,8]]*3,v_name=['V*01','V*01',None],j_name=['J*01']*3,c_name=['C*01',None,'C*01']))
    p=pl.DataFrame(dict(peptide_id=['p'],seq=['CA']))
    a=pl.DataFrame(dict(clone_id=['a','a','b'],peptide_id=['p']*3,start=[1]*3,llr=[1.,1.,-1.]))
    return SimpleNamespace(clonotypes=c,peptides=p,alignments=a),pl.DataFrame(dict(clone_id=['a','b','ref'],confirmed_nats=[20.,0.,0.]))


def test_union_partial_shared_junction_and_fixed_views():
    ev,s=fixture();out=capture_tables(ev,s,{'ref'},chunk_size=1)
    j=out['junctions'].filter((pl.col('view')=='all_compatible_alignments')&(pl.col('population')=='all'))
    assert j['covered_aa'].item()==2 and j['member_records'].item()==3
    assert j['any_partial_capture'].item() and not j['complete_capture'].item()
    r=out['records'].filter((pl.col('view')=='positive_llr_alignments')&(pl.col('region')=='anchor_inclusive_junction'))
    assert r['covered_aa'].to_list()==[2,0,0]
    proxy=out['genes'].filter((pl.col('view')=='all_compatible_alignments')&(pl.col('segment')=='V'))
    assert proxy['direct_region_records'].sum()==0  # label5 is not direct V support
    for name,frame in capture_tables(ev,s,{'ref'},chunk_size=100).items():
        assert frame.sort(frame.columns).equals(out[name].sort(frame.columns))


def test_real_counts_strict_background_only_and_reference_unmeasured():
    ev,s=fixture();counts=pl.DataFrame(dict(clone_id=['a','b'],count=[3,7]))
    out=capture_tables(ev,s,{'ref'},counts=counts,count_unit='reads',count_provenance='measured source')
    assert out['records'].filter(pl.col('is_reference'))['airr_count'].is_null().all()
    g=out['genes'].filter((pl.col('view')=='positive_llr_alignments')&(pl.col('population')=='background')&(pl.col('segment')=='V'))
    assert g['airr_capture_fraction'].item()==.3
    with pytest.raises(ValueError):capture_tables(ev,s,{'ref'},counts=counts.head(1),count_unit='reads',count_provenance='x')
    with pytest.raises(ValueError):capture_tables(ev,s,{'ref'},counts=counts,count_unit='protein',count_provenance='x')


def test_integer_precision_overflow_and_empty_reference_population():
    ev,s=fixture();out=capture_tables(ev,s,set())
    assert out['junctions'].filter(pl.col('population')=='reference').is_empty()
    empty=out['record_capture'].filter(pl.col('population')=='reference')
    assert empty['records'].sum()==0 and empty['captured_defined_record_fraction'].is_null().all()
    counts=pl.DataFrame({'clone_id':['a','b'],'count':[2**53+1,7]},schema_overrides={'count':pl.UInt64})
    out=capture_tables(ev,s,{'ref'},counts=counts,count_unit='reads',count_provenance='fixture')
    assert out['records'].filter(pl.col('clone_id')=='a')['airr_count'][0]==2**53+1
    assert out['records']['airr_count'].dtype==pl.UInt64
    for bad in [pl.Series('count',[1.,2.]),pl.Series('count',[2**63,1],dtype=pl.UInt64)]:
        with pytest.raises(ValueError):capture_tables(ev,s,{'ref'},counts=counts.with_columns(bad),count_unit='reads',count_provenance='fixture')
    c=out['record_capture'].filter((pl.col('population')=='background')&(pl.col('region')=='C')&(pl.col('view')=='all_compatible_alignments'))
    assert c['captured_defined_record_fraction'].item()==0


def test_native_protection_is_not_ndn_and_cdr3_uses_coordinates():
    ev,s=fixture()
    ev.clonotypes=ev.clonotypes.with_columns(pl.Series('region_ids',[
        [0,4,6,6,7,8],[0,4,6,6,7,8],[0,5,6,6,7,8]]))
    baseline=capture_tables(ev,s,{'ref'})
    ev.clonotypes=ev.clonotypes.with_columns(
        pl.Series('junction_core_status',['inferred_vj_excluded_core',
            'conservative_whole_junction_unresolved',None]),
        pl.Series('junction_core_start_nt',[3,None,None],dtype=pl.Int64),
        pl.Series('junction_core_end_nt',[9,None,None],dtype=pl.Int64),
        pl.Series('junction_core_protection_provenance',['inferred protection; not N-only',
            'unresolved protection; not N-only',None]))
    out=capture_tables(ev,s,{'ref'})
    native=out['records'].filter((pl.col('clone_id')!='ref')&
        (pl.col('view')=='all_compatible_alignments'))
    assert native.filter(pl.col('region')=='NDN').is_empty()
    cdr3=native.filter(pl.col('region')=='CDR3')
    assert cdr3['defined_aa'].to_list()==[2,2] and cdr3['covered_aa'].to_list()==[1,1]
    core=native.filter(pl.col('region')=='junction_protection')
    assert core['defined_aa'].to_list()==[2,2] and core['covered_aa'].to_list()==[1,1]
    assert core['junction_core_status'].to_list()==[
        'inferred_vj_excluded_core','conservative_whole_junction_unresolved']
    assert core['junction_core_start_nt'].to_list()==[3,None]
    assert core['junction_core_end_nt'].to_list()==[9,None]
    assert core['junction_core_protection_provenance'].null_count()==0
    # Supplied-query rows retain every legacy column/value; classes stay based
    # on explicit anchored junction sequence, independently of protection.
    assert out['records'].filter(pl.col('clone_id')=='ref').select(baseline['records'].columns).equals(
        baseline['records'].filter(pl.col('clone_id')=='ref'))
    assert out['junctions'].sort(out['junctions'].columns).equals(baseline['junctions'].sort(baseline['junctions'].columns))
    summary=out['record_capture'].filter((pl.col('population')=='background')&
        (pl.col('view')=='all_compatible_alignments')&(pl.col('region')=='junction_protection'))
    assert summary['records'].item()==2 and summary['defined_aa'].item()==4
    assert out['record_capture'].filter((pl.col('population')=='background')&
        (pl.col('region')=='NDN')).is_empty()


def test_null_native_status_keeps_all_legacy_outputs_unchanged():
    ev,s=fixture();baseline=capture_tables(ev,s,{'ref'})
    ev.clonotypes=ev.clonotypes.with_columns(junction_core_status=pl.lit(None,dtype=pl.String))
    for name,frame in capture_tables(ev,s,{'ref'}).items():
        assert frame.sort(frame.columns).equals(baseline[name].sort(baseline[name].columns))
    ev.clonotypes=ev.clonotypes.with_columns(junction_core_status=pl.lit('inferred_vj_excluded_core'))
    with pytest.raises(ValueError,match='bounds and provenance'):
        capture_tables(ev,s,{'ref'})


def test_explicit_unit_placeholders_cover_all_candidates_without_claiming_reads():
    ev,s=fixture();counts=ev.clonotypes.select('clone_id').with_columns(count=pl.lit(1,dtype=pl.UInt64))
    with pytest.raises(ValueError):
        capture_tables(ev,s,{'ref'},count_unit='placeholder',count_provenance='author-requested count=1')
    out=capture_tables(ev,s,{'ref'},counts=counts,count_unit='placeholder',count_provenance='author-requested count=1')
    r=out['records']
    assert r['airr_count'].min()==r['airr_count'].max()==1
    assert r['counts_are_placeholder'].all() and set(r['count_unit'])=={'placeholder'}
    total=out['record_capture'].filter((pl.col('view')=='all_compatible_alignments')&(pl.col('region')=='protein')&(pl.col('population')=='all'))
    assert total['airr_count'].item()==3 and total['captured_airr_count'].item()==2
    assert total['weighted_mean_coverage'].item()==pytest.approx(2/9)
    with pytest.raises(ValueError):
        capture_tables(ev,s,{'ref'},counts=counts.with_columns(count=pl.lit(2)),count_unit='placeholder',count_provenance='author-requested count=1')


def test_peptide_support_counts_species_rows_exact_and_zero_hit_candidates():
    from vdjaqp.report import peptide_support_table
    ev,_=fixture()
    ev.peptides=pl.DataFrame({'peptide_id':['p','q','r'],'seq':['IL','LL','AA']})
    ev.matches=pl.DataFrame({'clone_id':['a','a','a','a','b'],'peptide_id':['p','p','q','r','r'],'n_subs':[0,0,0,1,0]})
    table=peptide_support_table(ev).sort('clone_id')
    assert table.select('compatible_peptide_rows','compatible_peptide_species','exact_peptide_species').rows()==[(3,2,1),(1,1,1),(0,0,0)]
