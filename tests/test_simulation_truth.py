"""Truth regions preserve full denominators and native junction anchors."""
import pytest
from vdjaqp.simulation_truth import regions_from_annotations


def fixture(locus='IGH'):
    chain=dict(chain_id='c',locus=locus,airr_truth={'junction_aa':'CDEFGW','v_call':'V','j_call':'J'},constant_nt_gene='C')
    annotation={'sequence_id':'c','rev_comp':'F','junction_aa':'CDEFGW'}
    aa='AAAAACDEFGWTT';chain['variable_aa']=aa;chain['protein_aa']=aa+'PP';chain['variable_nt']='GCT'*5+'TGTGATGAATTTGGTTGGACTACT'
    chain['airr_truth'].update(junction=chain['variable_nt'][15:33],junction_start=15,junction_end=33)
    annotation['sequence']=chain['variable_nt']
    for name,start,end in [('fwr1',0,1),('cdr1',1,2),('fwr2',2,3),('cdr2',3,4),('fwr3',4,6),('cdr3',6,10),('fwr4',10,13)]:
        annotation[name+'_start']=start*3+1;annotation[name+'_end']=end*3
    return chain,annotation


@pytest.mark.parametrize('locus',['IGH','IGK','IGL'])
def test_partition_and_constants(locus):
    c,a=fixture(locus);metadata,regions=regions_from_annotations([c],[a])
    r={v['region']:v for v in regions.to_dicts()}
    assert metadata['locus'][0]==locus
    assert (r['JUNCTION']['start'],r['JUNCTION']['end'])==(5,11)
    assert r['FR3']['end']==6 and r['FR4']['start']==10
    assert (r['C']['start'],r['C']['end'])==(13,15)
    assert sum(v['end']-v['start'] for k,v in r.items() if k not in ('JUNCTION','GENERATOR_CDR3','V','D','J'))==15


def test_missing_and_invalid_boundaries():
    c,a=fixture();a['fwr1_start']=''
    _,r=regions_from_annotations([c],[a])
    assert r.filter(r['region']=='FR1')['status'][0]=='unavailable'
    assert r.filter(r['region']=='FR1')['start'][0] is None
    c,a=fixture();a['fwr1_end']=4
    _,r=regions_from_annotations([c],[a])
    assert r.filter(r['region']=='FR1')['status'][0]=='unavailable_non_codon_boundary'
    c,a=fixture();a['fwr1_end']=999
    with pytest.raises(ValueError):regions_from_annotations([c],[a])


def test_terminal_partial_codon_is_not_a_variable_residue():
    c,a=fixture()
    c['variable_nt'] += 'G'
    a['sequence']=c['variable_nt']
    a['fwr4_end']=len(c['variable_nt'])
    metadata,regions=regions_from_annotations([c],[a])
    assert metadata['untranslated_terminal_nt'][0]=='G'
    assert regions.filter(regions['region']=='FR4')['end'][0]==len(c['variable_aa'])
    a['fwr4_end']-=1
    regions_from_annotations([c],[a])
    a['fwr4_end']+=2
    with pytest.raises(ValueError):regions_from_annotations([c],[a])


def test_incomplete_annotation_keeps_truth_junction_without_filling_regions():
    c,a=fixture()
    a['junction_aa']=''
    a['fwr4_start']=a['fwr4_end']=''
    a['cdr3_end']=27
    _,r=regions_from_annotations([c],[a])
    by={v['region']:v for v in r.to_dicts()}
    assert by['CDR3']['status']=='unavailable_junction_boundary_conflict'
    assert by['CDR3']['start'] is None
    assert (by['JUNCTION']['start'],by['JUNCTION']['end'])==(5,11)
    c,a=fixture();c['protein_aa']='BAD'
    with pytest.raises(ValueError):regions_from_annotations([c],[a])


def test_empty_and_full_rna_validation():
    metadata,regions=regions_from_annotations([],[])
    assert metadata.height==regions.height==0
    c,a=fixture();c['rnaseq_nt']='TAA'
    with pytest.raises(ValueError):regions_from_annotations([c],[a])


def test_templated_segment_projection_and_gene_alias():
    c,a=fixture();c['airr_truth'].update(v_sequence_start=0,v_sequence_end=17,d_sequence_start=19,d_sequence_end=24,j_sequence_start=26,j_sequence_end=39,v_call='IGHV-A*01',d_call='IGHD-A*01',j_call='IGHJ2/3-A*01')
    _,regions=regions_from_annotations([c],[a]);r={v['region']:v for v in regions.to_dicts()}
    assert (r['V']['start'],r['V']['end'])==(0,5)
    assert (r['D']['start'],r['D']['end'])==(7,8)
    assert (r['J']['start'],r['J']['end'])==(9,13)
    assert r['J']['gene_call']=='IGHJ2/3-A*01'
    assert r['C']['gene_call']=='C' and r['FR1']['gene_call'] is None
    c,a=fixture('IGL');_,regions=regions_from_annotations([c],[a]);d=regions.filter(regions['region']=='D').row(0,named=True)
    assert d['start'] is None and d['end'] is None and d['status']=='not_applicable'


def test_reference_gene_label_is_additional_scalar_metadata():
    c,a=fixture();c.update(canonical_v_gene='IGHV1-18',canonical_v_match_status='unique_gene',canonical_v_gene_count=1)
    metadata,_=regions_from_annotations([c],[a])
    assert metadata['v_call'][0]=='V'
    assert metadata['canonical_v_gene'][0]=='IGHV1-18'
    assert metadata['canonical_v_gene_count'][0]==1
    assert 'IGHV1-18' in metadata.write_csv(separator='\t')
    c.update(canonical_v_gene=None,canonical_v_match_status='ambiguous_genes',canonical_v_gene_count=2)
    metadata,_=regions_from_annotations([c],[a])
    assert metadata['canonical_v_gene'][0] is None
    assert metadata['canonical_v_match_status'][0]=='ambiguous_genes'
    c['canonical_v_gene']='IGHV1-18'
    with pytest.raises(ValueError,match='canonical V'):
        regions_from_annotations([c],[a])


def test_component_latent_parent_and_shm_are_flat_evaluation_metadata():
    c,a=fixture()
    c.update(component='memory',ancestor_id='latent0',parent_id='latent0',
             branch_distance={'substitutions':2,'denominator_nt':39,'fraction':2/39},
             root_distance={'substitutions':3,'denominator_nt':39,'fraction':3/39})
    metadata,regions=regions_from_annotations([c],[a])
    row=metadata.row(0,named=True)
    assert row['component']=='memory' and row['ancestor_id']==row['parent_id']=='latent0'
    assert row['branch_substitutions']==2 and row['ancestor_substitutions']==3
    assert row['branch_denominator_nt']==row['ancestor_denominator_nt']==39
    assert row['ancestor_fraction']==3/39
    assert metadata.height==1 and regions['chain_id'].unique().to_list()==['c']
    assert 'latent0' in metadata.write_csv(separator='\t')


def test_internal_w_alternative_annotation_preserves_known_junction_and_cdr3():
    c,a=fixture()
    # Same physical native genotype, with an internal W before the true terminal W.
    c['variable_aa']='AAAAACDEWGWTT'
    c['protein_aa']=c['variable_aa']+'PP'
    c['variable_nt']='GCT'*5+'TGTGATGAATGGGGTTGGACTACT'
    c['airr_truth'].update(junction_aa='CDEWGW',junction=c['variable_nt'][15:33])
    a.update(sequence=c['variable_nt'],junction_aa='CDEW',junction=c['variable_nt'][15:27],
             cdr3_end=24,fwr4_start=25,fwr4_end=27,j_call='alternative_internal_J')
    metadata,regions=regions_from_annotations([c],[a])
    row=metadata.row(0,named=True)
    assert row['junction_annotation_concordance']=='conflict'
    assert row['annotated_junction_aa']=='CDEW' and row['generator_junction_aa']=='CDEWGW'
    assert row['annotated_j_call']=='alternative_internal_J' and row['j_call']=='J'
    assert row['generator_junction_nt']==c['airr_truth']['junction']
    by={r['region']:r for r in regions.to_dicts()}
    assert (by['JUNCTION']['start'],by['JUNCTION']['end'])==(5,11)
    assert (by['GENERATOR_CDR3']['start'],by['GENERATOR_CDR3']['end'])==(6,10)
    assert by['GENERATOR_CDR3']['source']=='generator_junction'
    for region in ('CDR3','FR4'):
        assert by[region]['start'] is None
        assert by[region]['status']=='unavailable_junction_boundary_conflict'
    assert metadata.height==1 and by['FR3']['end']==6


def test_native_junction_nt_corruption_or_coordinate_change_is_rejected():
    for changes in ({'junction_start':16}, {'junction_end':32},
                    {'junction':'TGTGATGAATTTGGTTGA'}, {'junction_start':None}):
        c,a=fixture()
        c['airr_truth'].update(changes)
        with pytest.raises(ValueError,match='native genotype junction'):
            regions_from_annotations([c],[a])
    c,a=fixture()
    meta,_=regions_from_annotations([c],[a])
    assert meta['junction_annotation_concordance'][0]=='matches'
    a['junction_aa']=''
    meta,_=regions_from_annotations([c],[a])
    assert meta['junction_annotation_concordance'][0]=='missing'


@pytest.mark.parametrize('partial_region',['FR1','FR4'])
def test_partial_inferred_terminal_regions_keep_all_known_truth(partial_region):
    c,a=fixture()
    if partial_region=='FR1':
        # A valid seven-region partition apart from the unaligned first codon.
        for name,start,end in [('fwr1',1,2),('cdr1',2,3),('fwr2',3,4),
                               ('cdr2',4,5),('fwr3',5,6)]:
            a[name+'_start']=start*3+1
            a[name+'_end']=end*3
    else:
        a['fwr4_end']=36
    metadata,regions=regions_from_annotations([c],[a])
    by={r['region']:r for r in regions.to_dicts()}
    assert by[partial_region]['start'] is None and by[partial_region]['end'] is None
    assert by[partial_region]['status']=='unavailable_incomplete_variable_partition'
    assert metadata['variable_aa'][0]==c['variable_aa']
    assert metadata['protein_aa'][0]==c['protein_aa']
    assert (by['JUNCTION']['start'],by['JUNCTION']['end'])==(5,11)
    assert (by['GENERATOR_CDR3']['start'],by['GENERATOR_CDR3']['end'])==(6,10)
    assert (by['C']['start'],by['C']['end'])==(13,15)
    assert all(by[r]['status']=='available' for r in ('JUNCTION','GENERATOR_CDR3','C'))
    # A partial outer annotation never licenses an internal partition gap.
    if partial_region=='FR4':
        a['fwr3_start']+=3
        with pytest.raises(ValueError,match='complete region partition'):
            regions_from_annotations([c],[a])
