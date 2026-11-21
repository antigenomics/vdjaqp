---
name: vdjaqp
description: Develop or run VDJ-AQP peptide-to-repertoire workflows, choose Evidence APIs, and verify placement and scoring contracts. Use within the VDJ-AQP repository or when preparing its candidate, peptide and model inputs.
---

# VDJ-AQP workflow

Read `AGENTS.md` and the relevant `docs/{pipeline,input-output,parameters}.rst`
from the repository root. These own current schemas and defaults.

- Supplied observations and an independently calibrated model:
  `Evidence.from_frames(candidates, peptides, error_model=model)` in
  `vdjaqp/pipeline.py`. `Evidence.build` additionally loads a named dataset and
  calibrates; do not substitute it for the explicit-model route.
- Mapping belongs to `PeptideMatcher` and seqtree batch search. Keep placements
  separate from credited evidence and accepted identities. Default search is
  ungapped; optional gapped spans/counts are search-only, rejected by residue scoring.
- Preserve I/L equivalence, zero-based half-open reference intervals, explicit
  `background_member`, measured sequence extent, missingness and ambiguous identity.
  CDR3 excludes conserved anchors; junction includes them.
- Resolve `vdjaqp.resolve_cpu_budget` before numerical imports. Native searches
  take batches within that ceiling; do not wrap them in per-query process pools.
- Nextflow and Snakemake observed analysis use `sample,candidates,peptides,producer,model` CSV inputs
  and fresh output paths. The producer is the actual upstream tool; adapters do
  not change provenance. Caching/resume are disabled.
- ARDA processes both repertoire geometries; MIGEC supplies supported UMI
  processing. Neither owns peptide mapping. See `docs/preparation.rst` for RNA
  expression and MS-raw provenance.

Verify software with the locked environment, pytest, the toy example and the
warning-as-error docs build in `CONTRIBUTING.md`. For a data-free Nextflow check,
create fresh inputs with `workflows/make_fixture.py`; its model is only a technical
fixture. Software verification does not authorize scientific runs, controls,
threshold changes or recalibration. For a paper rerun, retain the exact software
commit, dataset revision, input hashes, model and study-specific rule. Apply
the 20-confirmed-nat rule and the extra junction filter as separate comparisons;
use the declared RNA-adjusted removal or clone rule for their respective studies.
See `docs/data.rst` for published versus current runtime versions. Bootstrap
companion data before readers;
`vdjaqp-data` publication is per-file, not an atomic snapshot.

Run Snakemake from the repository root with a fresh output directory. CI compares
all observed-analysis Parquet tables across both workflow engines.

```sh
uvx --from snakemake==9.13.7 snakemake --cores 4 --config samples=samples.csv outdir=results/new-snakemake-analysis
```

Choose ARDA RNA-seq geometry for fragmented RNA reads and amplicon geometry
for AIRR libraries. For supported UMI-tagged AIRR libraries, assemble molecule
consensus before annotation. See `examples/inputs` for the preparation
worksheet, public tagged FASTQ and bundled mock peptide calls.
