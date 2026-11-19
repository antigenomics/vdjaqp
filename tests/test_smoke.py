"""Synthetic migration checks; no private datasets or scientific reruns."""
import polars as pl
from vdjaqp.data import to_gene_level, isotype_of
from vdjaqp.match import PeptideMatcher


def test_gene_normalization_and_isotype_are_label_independent():
    frame = pl.DataFrame({'v_name': ['IGHV1-2*02'], 'j_name': ['IGHJ4*01'], 'c_name': ['IGHG1*03']})
    normalized = to_gene_level(frame)
    assert normalized['v_name'].to_list() == ['IGHV1-2']
    assert normalized['j_name'].to_list() == ['IGHJ4']
    assert normalized.select(isotype_of(normalized)).to_series().to_list() == ['IgG']


def test_batched_search_preserves_mass_alphabet_and_empty_schema():
    matcher = PeptideMatcher(['clone0'], ['ACDEFGHIKLMNPQRSTVWY'])
    hits = matcher.search(['pep0', 'short'], ['GHIKLM', 'ACD'], max_subs=0, threads=1)
    assert hits.select('clone_id', 'peptide_id', 'start', 'end', 'n_subs').rows() == [('clone0', 'pep0', 5, 11, 0)]
    assert matcher.search([], [], threads=1).schema == hits.schema


def test_search_honors_process_budget_in_one_batch():
    import subprocess
    import sys
    code = """
import vdjaqp
vdjaqp.resolve_cpu_budget(1)
from vdjaqp.match import PeptideMatcher
class Index:
    calls = []
    def search_batch(self, queries, **kwargs):
        self.calls.append((queries, kwargs['threads']))
        return empty_result
matcher = PeptideMatcher(['clone0'], ['ACDEFGHIKLMNPQRSTVWY'])
empty_result = matcher._index.search_batch(['WWWWWW', 'WWWWWW'], max_subs=0, threads=1)
index = Index()
matcher._index = index
matcher.search(['p0','p1'], ['ACDEFG','HIKLMN'])
assert len(index.calls) == 1 and index.calls[0][1] == 1
"""
    subprocess.run([sys.executable, '-c', code], check=True, timeout=30)
