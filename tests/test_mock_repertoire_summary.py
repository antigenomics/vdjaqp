"""Constructed native metadata only; no scientific cohort or native tools."""
from copy import deepcopy
import gzip
import json
import polars as pl
import pytest
from scripts.summarize_mock_repertoire import summarize,main


def fixture():
    rows=[]
    for i,light in enumerate(['IGK','IGL']):
        for locus in ['IGH',light]:
            rows.append(dict(chain_id=f'chain{i}_{locus}',antibody_id=f'pair{i}',locus=locus,
                family_id='f0' if i==0 else None,component='related' if i==0 else None,
                variable_nt='GCTTGTGCTTTT',coding_variable_nt='GCTTGTGCTTTT',variable_aa='ACAF',protein_aa='ACAFG',
                constant_nt_gene={'IGH':'IGHG1','IGK':'IGKC','IGL':'IGLC2'}[locus],
                root_distance={'substitutions':i,'denominator_nt':12} if i==0 else None,
                branch_distance={'substitutions':0,'denominator_nt':12},
                airr_truth=dict(productive=True,junction='TGTGCTTTT',junction_aa='CAF',
                    v_call=locus+'V1*01',j_call=locus+'J1*01',n_mutations=0 if i==0 else None)))
    return {'chains':rows,'latent_ancestors':[{'chain_id':'ancestor0','observed':False}]}


def test_all_loci_gene_length_identity_and_missing_shm_are_explicit():
    out=summarize(fixture(),expected_chains=4)
    assert set(out['gene_usage']['locus'])=={'IGH','IGK','IGL'}
    for namespace in ['native_GenAIRR_allele_calls','native_GenAIRR_gene_calls']:
        part=out['gene_usage'].filter(pl.col('gene_namespace')==namespace)
        assert part['numerator'].sum()==12
        assert part.filter(pl.col('locus')=='IGH')['denominator'].to_list()==[2]*3
    genes=out['gene_usage'].filter((pl.col('gene_namespace')=='native_GenAIRR_gene_calls')&(pl.col('segment')=='V'))
    assert all('*' not in s for s in genes['gene_call'])
    assert out['identity_counts'].filter(pl.col('locus')=='IGH')['chain_records'].item()==2
    assert out['identity_counts'].filter(pl.col('locus')=='IGH')['anchored_junction_AA_classes'].item()==1
    assert out['family_sizes']['family_count'].sum()==2
    assert out['population_counts']['missing_family_records'].sum()==2
    unknown=out['shm_descriptors'].filter((pl.col('endpoint')=='root')&(pl.col('component')=='missing'))
    assert unknown['available_endpoint_records'].sum()==0 and unknown['substitutions'].null_count()==2
    assert unknown['denominator_nt'].null_count()==2
    h=out['shm_distributions'].filter((pl.col('quantity')=='native_n_mutations')&pl.col('missing_value'))
    assert h['count'].sum()==2
    assert out['length_distributions']['count'].sum()==24


@pytest.mark.parametrize('problem',['count','duplicate','latent','pair','distance','productive'])
def test_strict_observed_membership_and_native_field_schema(problem):
    truth=deepcopy(fixture())
    if problem=='count':truth['chains'].pop()
    if problem=='duplicate':truth['chains'][1]['chain_id']=truth['chains'][0]['chain_id']
    if problem=='latent':truth['latent_ancestors'][0]['chain_id']=truth['chains'][0]['chain_id']
    if problem=='pair':truth['chains'][1]['antibody_id']='unpaired'
    if problem=='distance':truth['chains'][0]['root_distance']['substitutions']=-1
    if problem=='productive':truth['chains'][0]['airr_truth']['productive']=False
    with pytest.raises(ValueError):summarize(truth,expected_chains=4)


def test_cli_requires_validated_full_cohort_and_fresh_destination(tmp_path):
    truth=tmp_path/'truth.json.gz'
    raw=json.dumps(fixture()).encode()
    with gzip.open(truth,'wb') as stream:stream.write(raw)
    receipt=tmp_path/'validation.json'
    receipt.write_text(json.dumps(dict(status='complete',observed_chains=4)))
    flags=['--truth-gzip',str(truth),'--validation-receipt',str(receipt),'--output',str(tmp_path/'new')]
    with pytest.raises(ValueError,match='full-cohort'):main(flags)
    assert not (tmp_path/'new').exists()
    (tmp_path/'new').mkdir()
    with pytest.raises(FileExistsError):main(flags)


def test_cli_rejects_changed_raw_identity_before_summary(tmp_path):
    truth=tmp_path/'truth.json.gz'
    with gzip.open(truth,'wb') as stream:stream.write(json.dumps(fixture()).encode())
    receipt=tmp_path/'validation.json'
    receipt.write_text(json.dumps(dict(status='complete',observed_chains=66080,
        truth=dict(raw_sha256='0'*64,raw_bytes=1,gzip_sha256='0'*64))))
    output=tmp_path/'new'
    with pytest.raises(ValueError,match='raw/gzip identity'):
        main(['--truth-gzip',str(truth),'--validation-receipt',str(receipt),'--output',str(output)])
    assert not output.exists()


def test_exact_junction_and_il_classes_keep_distinct_definitions():
    truth=fixture()
    for row,aa,nt in [(truth['chains'][0],'CIF','TGTATTTTT'),
                      (truth['chains'][2],'CLF','TGTCTTTTT')]:
        row.update(variable_nt='GCT'+nt,coding_variable_nt='GCT'+nt,
            variable_aa='A'+aa,protein_aa='A'+aa+'G')
        row['airr_truth'].update(junction=nt,junction_aa=aa)
    heavy=summarize(truth,expected_chains=4)['identity_counts'].filter(pl.col('locus')=='IGH')
    assert heavy['anchored_junction_AA_classes'].item()==2
    assert heavy['anchored_junction_NT_classes'].item()==2
    assert heavy['anchored_junction_AA_IL_classes'].item()==1
    assert heavy['variable_AA_IL_classes'].item()==1
