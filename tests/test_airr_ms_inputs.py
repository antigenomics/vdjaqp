"""Sequence-preserving private protein inputs require explicit physical choices."""
import json
import polars as pl
import pytest
from test_adapter import example
from vdjaqp.airr_ms_inputs import prepare_airr_ms_inputs


@pytest.mark.parametrize('physical', [['clone_1'],['clone_1','other']])
def test_background_source_role_is_independent_of_physical_membership(tmp_path,physical):
    from vdjaqp.airr_ms_inputs import attach_background_membership
    frame=pl.concat([example(),example().with_columns(clone_id=pl.lit('other'))])
    tagged=attach_background_membership(frame,inserted_reference_ids=['clone_1'])
    assert tagged.drop('background_member').equals(frame)
    paths=prepare_airr_ms_inputs(tagged,physical_clone_ids=physical,equal_protein_amounts=True,
        source_format='software fixture',output=tmp_path/'bundle')
    assert pl.read_parquet(paths['candidates']).equals(tagged)
    mapping=pl.read_csv(paths['id_mapping'],separator='\t')
    assert mapping['background_member'].to_list()==[False,True]
    assert mapping['physical'].to_list()==[True,'other' in physical]
    assert json.loads(paths['manifest'].read_text())['background_population']['excluded']==1
    for ids in [['missing'],['clone_1','clone_1']]:
        with pytest.raises(ValueError):attach_background_membership(frame,inserted_reference_ids=ids)
    with pytest.raises(ValueError,match='inserted references'):
        attach_background_membership(frame.with_columns(background_member=pl.lit(True)),inserted_reference_ids=['clone_1'])


def test_preserves_sequences_candidates_and_explicit_membership(tmp_path):
    frame = pl.concat([example(), example().with_columns(clone_id=pl.lit('other'))])
    paths = prepare_airr_ms_inputs(frame, physical_clone_ids=['other'],
        equal_protein_amounts=True, source_format='synthetic annotated export', output=tmp_path/'bundle')
    assert pl.read_parquet(paths['candidates']).equals(frame)
    assert paths['search_fasta'].read_text() == '>airr_00000001\nAACAWWWGG\n>airr_00000002\nAACAWWWGG\n'
    assert paths['proteins_fasta'].read_text() == '>airr_00000002 [# intensity=10000 #]\nAACAWWWGG\n'
    mapping = pl.read_csv(paths['id_mapping'], separator='\t')
    assert mapping['clone_id'].to_list() == ['clone_1','other']
    assert mapping['physical'].to_list() == [False,True]
    manifest = json.loads(paths['manifest'].read_text())
    assert manifest['physical_rows'] == 1 and manifest['candidate_rows'] == 2
    assert manifest['native_protein_intensity'] == manifest['total_native_input_intensity'] == 10000
    assert mapping['protein_relative_amount'].to_list() == [None, 1.]
    assert mapping['native_protein_intensity'].to_list() == [None, 10000.]
    assert manifest['input_snapshot_sha256'] == manifest['artifact_sha256']['candidates']
    before = paths['candidates'].read_bytes()
    with pytest.raises(FileExistsError):
        prepare_airr_ms_inputs(frame, physical_clone_ids=['other'], equal_protein_amounts=True,
            source_format='synthetic annotated export', output=tmp_path/'bundle')
    assert paths['candidates'].read_bytes() == before


@pytest.mark.parametrize('selected,accepted,sequence', [
    ([],True,'AACAWWWGG'), (['missing'],True,'AACAWWWGG'),
    (['clone_1','clone_1'],True,'AACAWWWGG'), ('clone_1',True,'AACAWWWGG'),
    (['clone_1'],False,'AACAWWWGG'), (['clone_1'],True,'AACAWWWG*'),
])
def test_rejects_invalid_inputs_before_writes(tmp_path, selected, accepted, sequence):
    frame = example().with_columns(seq=pl.lit(sequence))
    with pytest.raises(ValueError):
        prepare_airr_ms_inputs(frame, physical_clone_ids=selected, equal_protein_amounts=accepted,
            source_format='synthetic annotated export', output=tmp_path/'bundle')
    assert not (tmp_path/'bundle').exists()


@pytest.mark.parametrize('scale', [0, -1, float('nan'), float('inf'), True, '10000'])
def test_invalid_native_scale_never_writes(tmp_path, scale):
    with pytest.raises(ValueError, match='native protein intensity'):
        prepare_airr_ms_inputs(example(), physical_clone_ids=['clone_1'],
            equal_protein_amounts=True, source_format='synthetic annotated export',
            output=tmp_path/'bundle', native_protein_intensity=scale)
    assert not (tmp_path/'bundle').exists()


def test_explicit_native_scale_and_total(tmp_path):
    frame = pl.concat([example(), example().with_columns(clone_id=pl.lit('other'))])
    paths = prepare_airr_ms_inputs(frame, physical_clone_ids=['other','clone_1'],
        equal_protein_amounts=True, source_format='synthetic annotated export',
        output=tmp_path/'bundle', native_protein_intensity=20000)
    assert paths['proteins_fasta'].read_text().count('[# intensity=20000 #]') == 2
    assert json.loads(paths['manifest'].read_text())['total_native_input_intensity'] == 40000
