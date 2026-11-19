"""Constructed observations only; no cohort regeneration or scientific runs."""
import hashlib
import json
import polars as pl
import pytest
from test_contigs import fixture
from scripts.export_observed_repertoire import export_observations, main


def inputs(tmp_path, geometry='amplicon'):
    raw, ann, support = fixture()
    name = geometry + '_c'
    raw[0]['sequence_id'] = ann[0]['sequence_id'] = name
    row = dict(support['s_c'], contig_id='c')
    root = tmp_path/'input'
    (root/'assembly/read_support').mkdir(parents=True)
    for path, rows in [('raw_support.json', raw), ('arda_contigs.json', ann),
                       (f'read_support/{geometry}_support.json', [row])]:
        (root/'assembly'/path).write_text(json.dumps(rows))
    return root, raw, ann, row


def test_export_preserves_native_observations_and_membership(tmp_path):
    from vdjaqp.contigs import adapt_contigs
    root, raw, ann, support = inputs(tmp_path)
    output = tmp_path/'export'
    result = export_observations(root, 'amplicon', output, cpus=1)
    expected = adapt_contigs(raw, ann, {'amplicon_c': support})
    for key, native in expected.items():
        exported = pl.read_parquet(output/(key+'.parquet'))
        assert exported.drop('background_member').equals(native)
        assert exported['background_member'].to_list() == [True]*native.height
    assert result['counts']['candidates'] == 1
    assert result['observed_locus_counts'] == {'IGH': 1}
    assert result['cpu_resolved'] == 1
    assert json.loads((output/'raw.json').read_text()) == raw
    for source in result['sources'].values():
        from pathlib import Path
        assert source['sha256'] == hashlib.sha256(Path(source['path']).read_bytes()).hexdigest()
    with pytest.raises(FileExistsError):
        export_observations(root, 'amplicon', output)


def test_support_mask_and_ids_validated_by_native_adapter(tmp_path):
    root, raw, ann, row = inputs(tmp_path)
    support_file = root/'assembly/read_support/amplicon_support.json'
    row['retained_primary_exact_match_mask'] = []
    support_file.write_text(json.dumps([row]))
    with pytest.raises(ValueError, match='support mask'):
        export_observations(root, 'amplicon', tmp_path/'badmask')
    row['retained_primary_exact_match_mask'] = [True]*len(row['sequence'])
    row['contig_id'] = 'wrong'
    support_file.write_text(json.dumps([row]))
    with pytest.raises(ValueError, match='keys differ'):
        export_observations(root, 'amplicon', tmp_path/'badids')
    assert not (tmp_path/'badmask').exists()
    assert not (tmp_path/'badids').exists()


def test_rna_boundary_search_is_one_public_context_batch(tmp_path, monkeypatch):
    import vdjaqp.transcript_boundary as boundary
    root, raw, ann, row = inputs(tmp_path, 'rnaseq')
    contexts = [dict(target='public_context', sequence='ACGT', mature_start=2)]
    context_file = tmp_path/'contexts.json'
    context_file.write_text(json.dumps(contexts))
    import functools
    import vdjaqp.contigs as contigs
    original_adapter = contigs.adapt_contigs
    adapter_calls = []
    @functools.wraps(original_adapter)
    def adapt(records, annotations, supports, **kwargs):
        adapter_calls.append((records, annotations, supports, kwargs))
        return original_adapter(records, annotations, supports, **kwargs)
    monkeypatch.setattr(contigs, 'adapt_contigs', adapt)
    projected = pl.DataFrame(dict(sequence_id=['rnaseq_c'], status=['unresolved'], mature_start=[None]))
    calls = []
    def search(records, annotations, supports, catalogue, *, threads, alignment_output):
        calls.append((records, annotations, supports, catalogue, threads))
        return projected
    monkeypatch.setattr(boundary, 'search_transcript_boundaries', search)
    receipt = export_observations(root, 'rnaseq', tmp_path/'export', cpus=1,
                                  boundary_contexts=context_file)
    assert len(calls) == 1
    assert calls[0] == (raw, ann, {'rnaseq_c': row}, contexts, 1)
    assert len(adapter_calls) == 1
    assert adapter_calls[0][:3] == (raw, ann, {'rnaseq_c': row})
    assert adapter_calls[0][3] == {'boundaries': projected}
    assert receipt['boundary_contexts_used']


@pytest.mark.parametrize('argument', ['--truth', '--generation', '--model', '--native-offset'])
def test_truth_generation_arguments_rejected(argument):
    with pytest.raises(SystemExit) as exc:
        main(['--root', 'unused', '--output', 'unused', '--geometry', 'amplicon', argument, 'unused'])
    assert exc.value.code == 2


def test_nonheavy_observations_retained_without_forced_candidates(tmp_path):
    root, raw, ann, row = inputs(tmp_path)
    ann[0]['locus'] = 'IGK'
    (root/'assembly/arda_contigs.json').write_text(json.dumps(ann))
    output = tmp_path/'light'
    receipt = export_observations(root, 'amplicon', output, cpus=1)
    assert receipt['observed_locus_counts'] == {'IGK': 1}
    assert receipt['counts'] == {'catalog': 1, 'candidates': 0, 'search_entries': 1}
    assert pl.read_parquet(output/'candidates.parquet')['background_member'].dtype == pl.Boolean
