# Native read-placement support

Uses existing minimap2 2.31-r1302 and isolated SAMtools/HTSlib 1.25.
Run `sh envs/read_support/setup.sh [positive_cpu_budget]` after baseline setup.
It verifies the official release archive, builds bundled HTSlib with the shared
CPU resolver and a 300-second compile limit, and refuses existing targets.
No global upgrade or source patch is used. Native compiler, zlib, bzip2 and
lzma are prerequisites; curses and cloud/libcurl features are disabled.
Binary identities in identity.json describe this host, not portable hashes.
Use `.tool-envs/samtools/bin/samtools`; the checker accepts --samtools overrides.

One native `minimap2 -ax sr --secondary=no` batch per scenario maps fresh paired
reads against raw native contigs. Separate sequential samtools sort and mpileup
stages consume that output. Preset sr enables two I/O threads, so mapper -t uses
budget minus two; samtools sort -@ uses budget minus one. Budget below three is
rejected rather than silently oversubscribed. No concurrent mapper/sort pipeline.

Mpileup uses -aa and a raw-contig FASTA; explicit excluded flags are UNMAP,
SECONDARY,QCFAIL,DUP,SUPPLEMENTARY. Other native defaults are unchanged: base
quality/BAQ minimum13, mapping quality minimum0, depth cap8000, BAQ and paired-overlap
handling enabled, anomalous pairs excluded. These are software-fixture policies,
not scientifically calibrated thresholds. No tuning scan occurred.

Exact called-base support counts only . and , tokens after removing read-start
MAPQ markers and indel payloads using the documented native grammar. N reference
positions never become supported merely by depth. A positive count means a
native-retained primary-compatible placement, not uniquely attributable molecule,
true antibody origin, calibrated quality or full haplotype correctness. Shared
sequences/primary tie choices can leave compatible contigs unsupported.

Official documentation checked against the installed CLI:
- https://raw.githubusercontent.com/samtools/samtools/1.25/doc/samtools-mpileup.1
- https://lh3.github.io/minimap2/minimap2.html

Raw consensus/native positive-weight masks remain separate; native pseudocounts
are not used in these read-placement masks. Sorted BAM and per-base support JSON
are retained; SAM and textual pileup are consumed scratch and removed.

## Version selection and native diagnostic

SAMtools 1.24 has an upstream mpileup reference read-ahead cache lifetime bug.
Official 1.25 fixes the use-after-free ([PR 2364](https://github.com/samtools/samtools/pull/2364),
[source](https://raw.githubusercontent.com/samtools/samtools/1.25/bam_plcmd.c)).
The isolated published release replaces 1.24 for new runs; the original global
binary remains untouched. A single all-contig native mpileup on the retained
software-fixture RNA BAM checked 5,110 reference positions against FASTA:
all matched exactly and zero reference NUL bytes occurred. This diagnostic
did not regenerate simulated data or rerun scientific inference.

The archive/binary/compiler identity and compact diagnostic outcome are pinned
in identity.json. The setup retains raw source/archive and compact build logs;
consumed compiler objects are cleaned after verification. Pre-build project
footprint was 2,049,933,312 bytes, comfortably below 10 GB. Platform-dependent
system libraries must be rechecked when recreating this build elsewhere.
