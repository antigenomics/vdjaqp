<p align="center">
<img alt="VDJ-AQP logo: skull with antibody antlers and peptide chains" src="assets/vdjaqp.png" width="320">
</p>

<h1 align="center">VDJ-guided Antibody Quantification from Proteomes</h1>
<p align="center">
<a href="https://github.com/antigenomics/vdjaqp/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/antigenomics/vdjaqp/actions/workflows/ci.yml/badge.svg?branch=master"></a>
<a href="https://docs.isalgo.dev/vdjaqp/"><img alt="docs" src="https://github.com/antigenomics/vdjaqp/actions/workflows/docs.yml/badge.svg?branch=master"></a>
<img alt="Python 3.11 or newer" src="https://img.shields.io/badge/python-3.11%2B-blue">
<a href="LICENSE"><img alt="License" src="https://img.shields.io/badge/license-Academic%20%26%20Non--Commercial-green"></a>
</p>
<p align="center"><a href="https://docs.isalgo.dev/vdjaqp/quickstart.html">Quickstart</a> · <a href="https://docs.isalgo.dev/vdjaqp/">Documentation</a> · <a href="https://docs.isalgo.dev/vdjaqp/input-output.html">Output reference</a> · <a href="https://docs.isalgo.dev/vdjaqp/methods.html">Methods</a></p>

**VDJ-AQP links LC-MS/MS antibody peptides to B-cell receptor repertoire candidates.**
It keeps compatible peptide placements, credited sequence evidence and accepted
identities separate. Shared observations compete across candidates, residue
evidence is bounded, and ambiguous identities remain explicit.

The software is a Python package, runnable scripts and a Nextflow workflow.
[seqtree](https://github.com/antigenomics/seqtree) supplies peptide matching;
[ARDA](https://github.com/antigenomics/arda) processes both amplicon and fragmented
RNA sequencing geometries. [MIGEC](https://github.com/antigenomics/migec) provides
supported UMI molecule processing. Heavy, kappa and lambda chains retain separate
region, gene and junction readouts.

## Install

```sh
git clone https://github.com/antigenomics/vdjaqp.git
cd vdjaqp
sh setup.sh
.venv/bin/python scripts/check_environment.py
```

Python >=3.11 is required. The pinned runtime is recreated with uv from `uv.lock`;
Nextflow execution additionally requires Java 17–25 and Nextflow 25.10.2.
Optional simulation and native-tool recipes are in [`envs/`](envs/).
The Python import is `vdjaqp`; this repository is the canonical software source.

## Quick start

Run the data-free tutorial:

```sh
.venv/bin/python examples/toy/run_toy.py
```

The toy demonstrates shared and distinguishing peptide evidence using invented
sequences. Follow the [tutorial](https://docs.isalgo.dev/vdjaqp/tutorial.html) to
interpret the candidate scores.

For an analysis, prepare candidate and peptide tables using the
[input/output contracts](https://docs.isalgo.dev/vdjaqp/input-output.html), and
copy [`workflows/samples.example.csv`](workflows/samples.example.csv):

```sh
nextflow run main.nf --samples samples.csv --outdir results/new-analysis
```

Output directories must be new. Caching and resume are disabled. The workflow
retains calibration and source identities; scientific defaults are documented
in the [parameter reference](https://docs.isalgo.dev/vdjaqp/parameters.html).

## Prepare your inputs

Use MIGEC for supported UMI extraction, correction and molecule consensus, and
ARDA for repertoire annotation in both sequencing geometries. Preserve observed
sequence extent, molecule/read provenance and chain identity. CDR3 excludes the
conserved anchors; junction includes them. Do not fill missing observed sequence
from germline.

Quantify matched RNA expression independently, preserving reference versions,
counts, TPM and their normalization denominators. Transcript frequency, UMI
support and fitted peptide mixture fractions are different quantities from
protein abundance.

Preserve MS vendor raw files, acquisition metadata, converter versions and
centroiding settings. Record the protease, modifications, search database,
target/decoy strategy and calibrated confidence. Retain run/spectrum identifiers
and distinguish database-search from de novo calls. The
[preparation guide](https://docs.isalgo.dev/vdjaqp/preparation.html) covers these
upstream steps and their provenance requirements.

## Data

Companion inputs, calibration and retained results are maintained at
[isalgo/igprotseq](https://huggingface.co/datasets/isalgo/igprotseq). Access rights
and redistribution terms are separate from the software license. Downloads use
the dataset manifest and SHA-256 verification:

```sh
uv sync --frozen --extra data
.venv/bin/vdjaqp-data --help
```

See the [data reference](https://docs.isalgo.dev/vdjaqp/data.html) for groups,
authentication and pinning a dataset revision. This repository contains no
private assay inputs or stakeholder provenance.

## Documentation and development

Documentation follows Diataxis: tutorials teach a complete example; how-to guides
cover setup and preparation; reference pages define parameters and tables;
explanation pages describe the evidence model and its assumptions.

```sh
uv sync --frozen --extra test --extra docs --extra data
.venv/bin/python -m pytest tests -q
make -C docs html
```

CI checks the locked installation, synthetic tests, toy example, bootstrap CLI
and documentation with warnings as errors. Native-tool integration tests require
the optional pinned runtimes; they do not run scientific datasets in CI.

`master` is stable; `dev` integrates development. Branch features from `dev` and
submit pull requests there; promote verified changes from `dev` to `master`.
The manuscript and its host-specific analysis records remain separately maintained.

## Citation and license

Cite this repository and the accompanying VDJ-AQP publication when available.

ISALGO LAB Academic and Non-Commercial Research License v1.0.
See [`LICENSE`](LICENSE). Third-party tools retain their own licenses.
