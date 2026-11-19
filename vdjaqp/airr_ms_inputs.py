"""Private, sequence-preserving protein inputs for AIRR-conditioned MS simulation."""
from pathlib import Path
import hashlib
import io
import json
import math
from numbers import Real

import polars as pl

from .adapter import to_airr_compatible


def attach_background_membership(candidates: pl.DataFrame, *, inserted_reference_ids) -> pl.DataFrame:
    """Carry source population role before renaming; physical membership is separate."""
    if isinstance(inserted_reference_ids,(str,bytes)):
        raise ValueError('inserted_reference_ids must be an explicit collection')
    references=list(inserted_reference_ids)
    if len(set(references))!=len(references) or not set(references)<=set(candidates['clone_id']):
        raise ValueError('inserted reference IDs must be unique and present')
    if 'background_member' in candidates.columns:
        role=candidates['background_member']
        if role.dtype!=pl.Boolean or role.null_count():
            raise ValueError('background_member must be Boolean and nonmissing')
        if candidates.filter(pl.col('clone_id').is_in(references))['background_member'].any():
            raise ValueError('inserted references cannot be background members')
        return candidates
    return candidates.with_columns((~pl.col('clone_id').is_in(references)).alias('background_member'))


def prepare_airr_ms_inputs(candidates: pl.DataFrame, *, physical_clone_ids,
                          equal_protein_amounts: bool, source_format: str,
                          output: Path, native_protein_intensity: float = 10000.,
                          protein_relative_amounts: dict | None = None) -> dict[str, Path]:
    """Write a fresh private input bundle; never infer physical membership.

    Equal amounts must be explicitly accepted. They are relative protein input
    weights, not RNA counts or measured serum abundances. All supplied sequence
    bytes, annotations, duplicate sequences and candidate rows are retained.
    """
    to_airr_compatible(candidates, source_format=source_format)
    if not candidates.height or candidates.select(
        (~pl.col('seq').str.contains('^[ACDEFGHIKLMNPQRSTVWY]+$')).fill_null(True).any()
    ).item():
        raise ValueError('all candidates must have nonempty canonical protein sequences')
    if protein_relative_amounts is None and equal_protein_amounts is not True:
        raise ValueError('equal protein amounts must be explicitly accepted')
    if protein_relative_amounts is not None and equal_protein_amounts is not False:
        raise ValueError('explicit amounts require equal_protein_amounts=False')
    if isinstance(physical_clone_ids, (str, bytes)):
        raise ValueError('physical_clone_ids must be an explicit nonempty sequence of IDs')
    selected = list(physical_clone_ids)
    if (not selected or any(not isinstance(x, str) for x in selected)
            or len(set(selected)) != len(selected)
            or not set(selected) <= set(candidates['clone_id'])):
        raise ValueError('physical IDs must be nonempty, unique and present in candidates')
    if (isinstance(native_protein_intensity, bool) or not isinstance(native_protein_intensity, Real)
            or not math.isfinite(native_protein_intensity) or native_protein_intensity <= 0
            or not math.isfinite(native_protein_intensity * len(selected))):
        raise ValueError('native protein intensity and total must be finite and positive')
    native_protein_intensity = float(native_protein_intensity)
    selected = set(selected)
    amounts = {name:1. for name in selected} if protein_relative_amounts is None else protein_relative_amounts
    if (not isinstance(amounts, dict) or set(amounts)!=selected
            or any(isinstance(v,bool) or not isinstance(v,Real) or not math.isfinite(v) or v<=0 for v in amounts.values())):
        raise ValueError('relative amounts require exactly the selected IDs and finite positive values')
    if protein_relative_amounts is not None and (len(selected)!=12 or not math.isclose(sum(amounts.values()),12.,rel_tol=0,abs_tol=1e-10)):
        raise ValueError('explicit panel amounts must contain twelve records and sum to twelve')
    intensities={name:float(v)*native_protein_intensity for name,v in amounts.items()}
    if any(not math.isfinite(v) or v<=0 for v in intensities.values()):
        raise ValueError('native per-record intensities must be finite and positive')
    mapping = candidates.select(['clone_id']+(['background_member'] if 'background_member' in candidates.columns else [])).with_columns(
        pl.Series('entry_id', [f'airr_{i+1:08d}' for i in range(candidates.height)]),
        pl.col('clone_id').is_in(selected).alias('physical'),
    ).with_columns(
        pl.Series('protein_relative_amount',[amounts.get(name) for name in candidates['clone_id']],dtype=pl.Float64),
        pl.Series('native_protein_intensity',[intensities.get(name) for name in candidates['clone_id']],dtype=pl.Float64))
    snapshot = io.BytesIO()
    candidates.write_parquet(snapshot, compression='zstd')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    paths = {key: output/name for key, name in {
        'candidates': 'candidates.parquet', 'search_fasta': 'search.fasta',
        'proteins_fasta': 'proteins.fasta', 'id_mapping': 'private_id_mapping.tsv',
        'manifest': 'manifest.json',
    }.items()}
    paths['candidates'].write_bytes(snapshot.getvalue())
    mapping.write_csv(paths['id_mapping'], separator='\t')
    with paths['search_fasta'].open('w') as search, paths['proteins_fasta'].open('w') as proteins:
        for entry, seq, physical, intensity in zip(mapping['entry_id'], candidates['seq'], mapping['physical'],mapping['native_protein_intensity']):
            search.write(f'>{entry}\n{seq}\n')
            if physical:
                proteins.write(f'>{entry} [# intensity={intensity:g} #]\n{seq}\n' if protein_relative_amounts is None else
                               f'>{entry} [# intensity={intensity:.17g} #]\n{seq}\n')
    manifest = dict(
        source_format=source_format, input_snapshot_sha256=hashlib.sha256(snapshot.getvalue()).hexdigest(),
        input_hash_definition='SHA256 of unchanged dataframe serialized as candidates.parquet; not original source-file bytes',
        candidate_rows=candidates.height, physical_rows=len(selected),
        sequence_policy='supplied amino-acid sequences unchanged, including any constant stub; no completion or reannotation',
        physical_membership='explicit caller-selected clone IDs; no inferred membership',
        protein_amount_assumption='equal relative protein amounts of 1.0 per selected record; not observed abundance or RNA counts',
        native_protein_intensity=native_protein_intensity,
        total_native_input_intensity=native_protein_intensity * len(selected),
        native_intensity_units='native MSSimulator input intensity; not experimentally calibrated; not normalized pool weights',
        native_default_source='https://github.com/OpenMS/OpenMS/blob/60ce8966d577b43a6419403dec171539e08fe513/src/utils/MSSimulator.cpp#L203-L230',
        confidentiality='all artifacts private; neutral IDs do not authorize release',
        artifact_sha256={key: hashlib.sha256(path.read_bytes()).hexdigest() for key,path in paths.items() if key != 'manifest'},
    )
    if protein_relative_amounts is not None:
        manifest.update(protein_amount_assumption='explicit relative panel weights summing to twelve; not measured abundance',
                        relative_amounts_sha256=hashlib.sha256(json.dumps(amounts,sort_keys=True).encode()).hexdigest(),
                        total_native_input_intensity=sum(intensities.values()))
    if 'background_member' in candidates.columns:
        manifest['background_population']=dict(source='explicit background_member candidate column; independent of physical membership',
            members=int(candidates['background_member'].sum()),excluded=candidates.height-int(candidates['background_member'].sum()))
    paths['manifest'].write_text(json.dumps(manifest, indent=2)+'\n')
    return paths
