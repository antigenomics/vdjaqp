# Input examples

[raw-inputs.tsv](raw-inputs.tsv) is a preparation worksheet. Replace the example
paths with files for the same specimen. Paired FASTQ records must be in matching
order. Leave read2 empty for single-end input and molecule_protocol empty for
libraries without UMIs. The protocol file records barcode coordinates and read
orientation from the library design.

| Library | Candidate preparation | Support unit |
| --- | --- | --- |
| Fragmented RNA-seq | ARDA RNA-seq geometry | Reads or original reported support |
| AIRR-seq without UMIs | ARDA amplicon geometry or annotated AIRR export | Reads |
| AIRR-seq with UMIs | MIGEC consensus for a supported protocol, then ARDA amplicon geometry | UMIs |

For each route, search the MS acquisition upstream and export peptide calls to
the VDJ-AQP Parquet schema. Score using the shared five-column sample sheet in
[workflows](../../workflows/samples.example.csv). A calibrated error model is
required for each analysis.

## Offline scoring fixture

The bundled gzip contains 10,000 invented database-search calls. Each row
represents a distinct mock spectrum. Generate the prepared Parquet inputs and
run the same workflow used with prepared user data:

```sh
.venv/bin/python workflows/make_fixture.py tmp/fixture --mock-ms
nextflow run main.nf --samples tmp/fixture/samples.csv --outdir tmp/fixture-nextflow
uvx --from snakemake==9.13.7 snakemake --cores 4 --config samples=tmp/fixture/samples.csv outdir=tmp/fixture-snakemake
```

The fixture exercises scoring and file formats. Its invented calls are not
instrument spectra or an empirical sensitivity benchmark.

## Public read fixture

[sc5p_v2_hs_PBMC_1k_t_cells1pct.fq.gz](sc5p_v2_hs_PBMC_1k_t_cells1pct.fq.gz)
contains 18,036 checked-out cDNA reads from 10x Genomics'
`sc5p_v2_hs_PBMC_1k` VDJ-T library. The upstream fixture selects whole cells
and retains cell-barcode and UMI tags. Use it to inspect the tagged FASTQ
format and test molecule processing. It is TCR data and has no matched
antibody MS acquisition.

The bytes match the [isalgo/umi_data fixture](https://huggingface.co/datasets/isalgo/umi_data/tree/23363ad9e7a248b834f175b31c3eb25349a0bf23/ci),
distributed upstream under [GPL-3.0](LICENSE.GPL-3.0). Upstream [source and selection recipe](https://huggingface.co/datasets/isalgo/umi_data/blob/23363ad9e7a248b834f175b31c3eb25349a0bf23/SOURCES.md).
SHA-256: `f14e8ec7fb251297eadcf4bf49840dcfb14d09b039d1b0de46149aabc8514f58`.

## AIRR Rearrangement example

[rearrangements.tsv](rearrangements.tsv) contains an invented AIRR record.
After making the offline MS fixture, append this line to
`tmp/fixture/samples.csv` and run either engine with a fresh output root:

```csv
airr,examples/inputs/rearrangements.tsv,tmp/fixture/peptides.parquet,synthetic AIRR,tmp/fixture/model.json
```

Use the repository root as the invocation directory. The candidate field
also accepts the same TSV compressed as `.tsv.gz`.
