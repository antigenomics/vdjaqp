# Published paired-read preprocessing

fastp 0.24.0, Bioconda osx-64 build `h9ea9c2a_1`, installed and executed through macOS Rosetta. Exact package URLs are in `osx-64-explicit.txt`; package hashes/license metadata and executable SHA256 are in `runtime.json`. This is a separate tool environment, not a change to the scientific scoring environment.

```sh
conda create --prefix .tool-envs/preprocess --file envs/preprocess/osx-64-explicit.txt -y
.tool-envs/preprocess/bin/fastp --version
```

Only conda-forge/Bioconda were used; no default-channel agreement or global configuration was changed. The native merge interface is `--merge --include_unmerged --merged_out FILE`: it writes merged fragments and both unmerged mates to one stream. Adapter trimming, poly-G trimming, quality filtering and length filtering are explicitly disabled for the generated adapter-free fixture. This is an overlap-merging stage, not an additional scientific filter or homemade sequencing-error model. The original amplicon/fragmentation scenario labels and ARDA modes remain distinct.

`smoke_preprocess.py` verifies native passed-read counts, native merged-read tags, output FASTQ record counts and the identity `output_records = merged_pairs + 2*unmerged_pairs`. A direct two-pair software fixture produced one merged fragment plus two unmerged mates, with all four input reads accounted for. Existing output is refused. ARDA then consumes this observed merged/unmerged stream; no true or germline sequence is supplied to its read reconstruction.

Publication: Chen et al., **fastp: an ultra-fast all-in-one FASTQ preprocessor.** Bioinformatics 2018, PMID [30423086](https://pubmed.ncbi.nlm.nih.gov/30423086/), DOI [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560). Primary documented merge semantics: [official fastp README](https://github.com/OpenGene/fastp#merge-paired-end-reads), additionally verified against installed `fastp --help`.

The previous ARDA 2.36 `--reconstruct` option cannot combine with its default junction-quality correction mode. The published fastp merger produces a genuine merged quality string and avoids inventing one in project code. This addition does not alter ARDA biological/search thresholds.

## Linux x86_64 artifact lock

`linux-64-explicit.txt` and `linux-64-artifacts.json` pin fastp 0.24.0 and minimap2 2.31; 13 artifacts, 42,896,956 advertised download bytes. Resolution used a metadata-only conda dry run with glibc 2.31. The JSON records SHA256, MD5, package sizes and dependency metadata; these sizes are archive downloads, not installed footprint. Linux x86_64/glibc 2.31 installation passed: all artifact archive SHA256/size checks, exact target versions, and native help/INI probes (OpenMS) or native version probes (fastp/minimap2). End-to-end Linux sequence/MS simulation remains unverified. Recreate a new Linux prefix with the platform-specific explicit lock; never copy this host's environment directory.
