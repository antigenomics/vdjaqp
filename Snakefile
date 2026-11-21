"""Run observed VDJ-AQP analysis from the Nextflow sample-sheet format."""

from pathlib import Path
import sys

script_root = str(Path(workflow.basedir) / "workflows")
sys.path.insert(0, script_root)
from sample_sheet import read_samples

samples = read_samples(config["samples"])
outdir = str(Path(config["outdir"]).resolve())
if any(character in outdir for character in "{}"):
    raise ValueError("Snakemake output paths cannot contain braces")
python = str(Path(workflow.basedir) / ".venv/bin/python")
driver = str(Path(workflow.basedir) / "scripts/workflow_analysis.py")


onstart:
    Path(outdir).mkdir(parents=True, exist_ok=False)


rule all:
    input:
        expand(outdir + "/{sample}/analysis", sample=list(samples))


rule analyse:
    input:
        candidates=lambda wc: samples[wc.sample]["candidates"],
        peptides=lambda wc: samples[wc.sample]["peptides"],
        model=lambda wc: samples[wc.sample]["model"],
    output:
        directory(outdir + "/{sample}/analysis"),
    params:
        producer=lambda wc: samples[wc.sample]["producer"],
        parent=lambda wc: str(Path(outdir) / wc.sample),
    wildcard_constraints:
        sample="[A-Za-z0-9_-]+",
    threads: 4
    resources:
        mem_mb=4096,
        runtime=15,
    shell:
        "mkdir -p -- {params.parent:q} && "
        "{python:q} {driver:q} {input.candidates:q} {input.peptides:q} "
        "{input.model:q} {output:q} --producer {params.producer:q} --cpus {threads}"
