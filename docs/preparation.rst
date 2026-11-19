Prepare repertoire, expression and MS inputs
============================================

Keep a private sample manifest connecting library, specimen, acquisition and
processing identifiers. Join modalities by a verified specimen key; where only
condition-level matching exists, report that level explicitly. Preserve raw files,
checksums, reference versions and complete commands separately from derived data.

Repertoire and molecule processing
----------------------------------

Use `MIGEC <https://github.com/antigenomics/migec>`_ for supported UMI molecule
processing and `ARDA <https://github.com/antigenomics/arda>`_ for both amplicon and
fragmented RNA sequencing geometries. Select the mode from the physical library
protocol; UMI positions, primers, read orientation and cell barcodes must come
from the assay metadata. Do not treat read counts as molecule counts. Assemble
libraries independently before constructing a qualified joint candidate union.
See :doc:`pipeline` and the pinned recipes in ``envs/`` for the implemented route.

Preserve observed variable/constant sequence and explicit missing extent; never
complete an observed sequence from germline. Keep heavy, kappa and lambda readouts
separate. CDR3 excludes conserved anchors; junction includes them. Annotation,
error correction and representative selection must retain their own provenance.

Expression from RNA-seq
-----------------------

Quantify matched whole-transcriptome RNA independently of peptide identification.
For bulk RNA-seq, `Salmon's documentation
<https://salmon.readthedocs.io/en/latest/salmon.html>`_ describes transcript
quantification and recommends a decoy-aware transcriptome. Record the transcript
FASTA and annotation release, index recipe, library strandedness, software version
and quantification command. Retain transcript identifiers and the mapping to genes;
keep counts and TPM as distinct quantities with named normalization denominators.

A reconstructed receptor sequence is not a measured transcript abundance. Generic
reference transcript estimates do not establish clonotype-specific expression.
UMI molecule support, read support, transcript abundance and protein abundance
remain different measurements. Do not infer serum concentration or binding titer
from RNA frequencies or fitted peptide mixture fractions. This is upstream
preparation guidance, not an expression quantifier implemented by VDJ-AQP.

MS raw data and peptide calls
-----------------------------

Preserve vendor raw files and acquisition metadata. Where an upstream engine
requires mzML, use a documented converter such as `ProteoWizard msconvert
<https://www.proteowizard.org/tools/msconvert.html>`_; record the converter version,
complete command, output checksum, centroiding and other filters. Choose those
settings for the acquisition and search engine rather than applying an arbitrary
universal conversion recipe. Retain the original raw files.

Record protease specificity, missed cleavages, modifications, mass tolerances,
instrument/acquisition mode, search engine/version, target/decoy construction and
search database checksum. Match protease metadata to the physical digest. Preserve
run and spectrum identifiers so alternative calls from one spectrum remain one
observation. Keep database-search and de novo calls and their confidence semantics
explicit; a search-engine q-value is not automatically the VDJ-AQP error model.

Use the measured peptide sequence, modifications/mass and calibrated confidence
fields required by :doc:`input-output`. Preserve I/L equivalence and distinguish
search compatibility, credited evidence and accepted identity. Do not use assay
truth, source labels or expected expression as scoring inputs. Raw conversion,
expression quantification and database search are upstream operations; the
Nextflow analysis entry consumes the documented candidate and peptide tables.

Data and reproducibility
------------------------

The companion dataset is `isalgo/igprotseq
<https://huggingface.co/datasets/isalgo/igprotseq>`_. See :doc:`data` for the
manifest-based bootstrap, access requirements and provenance. Dataset access does
not authorize redistribution of private inputs. Use fresh output directories;
retain software, input and model identities with each authorized analysis.
