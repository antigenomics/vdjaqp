# Pinned ARDA source overlay

The historical simulation overlay uses an isolated `arda-mapper==2.36.0` wheel.
The main scoring environment uses the version pinned in `uv.lock`. The overlay retains the
four-module AIRR projection/gzip-statistics backport from upstream source commit
`dc97fcfeb68e223a1688428e0fad38394449c3e8`, based on tag2.36.0 commit
`049553f35114f2578c99355f8afd4c107c42c4d9`. Its original focused21checks passed
before the subsequent native-membership addition; this does not establish
full-cohort feasibility. Historical outputs remain labelled by their producers.

`arda-2.36.0-read-membership.patch` additionally passes the EXISTING native
`correct_airr(..., read_map=...)` API through the shared pipeline exactly once,
writing `sequence_id -> clone_row` in `<prefix>.read_ids.tsv`. It does not change
correction/assembly defaults or insert any bases. The exporter uses native
corrected-clone row membership, preserving observed sequence/junction/V/J/C
calls separately from corrected-root calls. This avoids dropping longer raw
observations merely because their junction nucleotide keys were corrected.
Every current exporter call requires this native map and fails if it is absent.
Archived outputs remain readable without re-exporting them through an obsolete route.

Set up the main environment, then run from the project root:

```
sh setup.sh
sh envs/arda/setup_overlay.sh
```

The executable setup materialises the previously documented source recipe,
verifies upstream commit/license/patches and all five baseline/patched module
hashes, and leaves the main environment unchanged. Native extensions and vendor
libraries link to the separately installed baseline wheel; rebuild per host, never
transfer this directory. Existing overlay destinations are rejected.
The new destination is `.tool-envs/arda-overlay-read-map/overlay`; the earlier
four-module overlay is preserved separately. Nextflow supplies this PYTHONPATH
automatically for native FASTQ branches. Population runs verify all five actual
imported module hashes and record them with the source-patch identities.

`source-manifest.json` is authoritative. Preserve GPLv3 (`LICENSE.arda`, exact
upstream license), attribution and both patch files. No upstream release or
public dataset publication is implied. No stakeholder data are included.

The constructed synonymous-member regression and fresh actual both-preset
qualification are recorded in `docs/pipeline.rst`; do not infer
full native validation from the prior projection tests. Current resource and
cleanup bounds follow that finite specification and docs/pipeline.rst. No cache/resume
path or unrequested scientific rerun is authorised by this installation recipe.
