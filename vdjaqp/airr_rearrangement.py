"""Read AIRR Rearrangement TSV or TSV.gz into observed candidate tables."""

from pathlib import Path
import re

import polars as pl
from arda.refbuild.translate import translate

from .adapter import _validate
from .data import UNANNOTATED_REGION


def _unique_offset(sequence, query):
    start = sequence.find(query)
    if not query or start < 0 or sequence.find(query, start + 1) >= 0:
        raise ValueError("sequence must contain the annotated interval exactly once")
    return start


def _ungap(value):
    return (value or "").upper().replace(".", "").replace("-", "")


def _protein(row, junction_nt):
    amino = _ungap(row.get("sequence_alignment_aa"))
    if amino:
        return amino
    nt = _ungap(row.get("sequence_alignment"))
    if not re.fullmatch("[ACGT]+", nt) or not junction_nt:
        raise ValueError("provide sequence_alignment_aa or a nucleotide alignment and in-frame junction")
    phase = _unique_offset(nt, junction_nt) % 3
    coding = nt[phase:]
    return translate(coding[:len(coding) // 3 * 3])


def _candidate(row):
    junction_nt = (row.get("junction") or "").upper()
    translated = None
    if re.fullmatch("[ACGT]+", junction_nt) and len(junction_nt) % 3 == 0:
        translated = translate(junction_nt)
    seq = _protein(row, junction_nt if translated is not None else "")
    junction = (row["junction_aa"] or translated or "").upper()
    if translated is not None and translated != junction:
        raise ValueError("junction and junction_aa disagree")
    if not re.fullmatch("[A-Z]+", seq) or not re.fullmatch("C[A-Z]+[FW]", junction):
        raise ValueError("require an observed stop-free protein and anchor-inclusive junction_aa")
    start = _unique_offset(seq, junction)
    end = start + len(junction)
    labels = [UNANNOTATED_REGION] * len(seq)
    labels[start + 1:end - 1] = [5] * (len(junction) - 2)
    cursor = 0
    for label, field in enumerate(("fwr1_aa", "cdr1_aa", "fwr2_aa", "cdr2_aa", "fwr3_aa")):
        region = _ungap(row.get(field))
        if region:
            at = cursor + _unique_offset(seq[cursor:start + 1], region)
            cursor = at + len(region)
            labels[at:cursor] = [label] * len(region)
    fwr4 = _ungap(row.get("fwr4_aa"))
    if fwr4:
        at = end - 1 + _unique_offset(seq[end - 1:], fwr4)
        labels[at:at + len(fwr4)] = [7] * len(fwr4)
    if row.get("cdr3_aa") and row["cdr3_aa"].upper() != junction[1:-1]:
        raise ValueError("cdr3_aa must match junction_aa without its anchors")
    membership = row.get("background_member", "T")
    if membership not in ("T", "F"):
        raise ValueError("background_member must be T or F")
    return {
        "clone_id": row["sequence_id"], "seq": seq, "cdr3_seq": junction,
        "cdr3_start": start, "cdr3_end": end, "region_ids": labels,
        "v_name": row["v_call"], "j_name": row["j_call"], "c_name": row.get("c_call"),
        "background_member": membership == "T",
    }


def read_rearrangements(path: str | Path) -> pl.DataFrame:
    """Project AIRR records to the scoring schema without filtering the input."""
    source = pl.read_csv(path, separator="\t", infer_schema=False, quote_char=None)
    required = {"sequence_id", "v_call", "j_call", "junction_aa"}
    missing = required - set(source.columns)
    if missing:
        raise ValueError(f"missing AIRR columns: {sorted(missing)}")
    if source.is_empty():
        raise ValueError("empty AIRR Rearrangement table")
    candidates = []
    for index, row in enumerate(source.iter_rows(named=True), 2):
        try:
            candidates.append(_candidate(row))
        except ValueError as error:
            raise ValueError(f"AIRR row {index} ({row['sequence_id']!r}): {error}") from error
    result = pl.DataFrame(candidates, schema_overrides={
        "clone_id": pl.String, "v_name": pl.String, "j_name": pl.String, "c_name": pl.String,
    })
    _validate(result)
    return result
