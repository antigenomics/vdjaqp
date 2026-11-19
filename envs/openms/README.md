# OpenMS 3.0 isolated tool runtime

MSSimulator was removed in OpenMS 3.1. This macOS-arm64 host can run the official Intel macOS package through Rosetta; Docker daemon is absent. No global tool installation is required.

```sh
conda create --prefix .tool-envs/openms --file envs/openms/osx-64-explicit.txt --yes
```

Initial resolution: `conda create --prefix .tool-envs/openms --platform osx-64 --override-channels -c conda-forge -c bioconda openms=3.0.0=haa56855_1 --yes`. The explicit export pins transitive artifacts and hashes. Linux deployments must use a separately resolved linux-64 lock; never reuse this platform's environment directory.

Binary reports 3.0.0-pre-exported-20230807. Installed package recipe/version, rather than the date string alone, identifies this official release build.

## Linux x86_64 artifact lock

`linux-64-explicit.txt` and `linux-64-artifacts.json` pin OpenMS 3.0.0; 122 artifacts, 326,405,133 advertised download bytes. Resolution used a metadata-only conda dry run with glibc 2.31. The JSON records SHA256, MD5, package sizes and dependency metadata; these sizes are archive downloads, not installed footprint. Linux x86_64/glibc 2.31 installation passed: all artifact archive SHA256/size checks, exact target versions, and native help/INI probes (OpenMS) or native version probes (fastp/minimap2). End-to-end Linux sequence/MS simulation remains unverified. Recreate a new Linux prefix with the platform-specific explicit lock; never copy this host's environment directory.

## Required bounded shared-peptide simulator

After installing the pinned runtime, build the existing metadata-only fix:

```
.venv/bin/python scripts/build_openms_digest_fix.py --prefix "$PWD/.tool-envs/openms"
```

The baseline `MSSimulator` stays unchanged. Default project adapters now require
`MSSimulator-digestfix` and its verified build manifest, preventing silent use of
the quadratic repeated-accession metadata path. Sources, ABI, compiler, linked
libraries and executable hashes are recorded. The native fixture verifies the
same peptides/intensities/unique-accession sets and fixes duplicate evidence
rows6→3. No digestion/intensity/random model is changed. An explicit executable
override remains recorded as an override. The reproduction bootstrap builds
this scoped executable automatically; no new upstream package is released.
