Pipeline branches
=================

Run from the repository root with Java 17–25 and Nextflow **25.10.2**. ``sh setup.sh`` recreates the minimal Python runtime from ``uv.lock`` (no optional simulator tools). A local Nextflow launcher can be downloaded from https://github.com/nextflow-io/nextflow/releases/download/v25.10.2/nextflow-25.10.2-dist and marked executable. No global Python installation is needed.

Observed analysis
-----------------

::

   nextflow run main.nf --samples samples.csv --outdir results/new-analysis

CSV columns: ``sample,candidates,peptides,producer,model``. Use absolute input paths; sample names must be unique alphanumeric/underscore/hyphen identifiers. Shell-boundary path arguments accept letters, digits, spaces and ``_.,/:+=-`` only. Apostrophes and other shell metacharacters fail before task submission. ``producer`` states the actual producer (e.g. MiXCR), never ARDA merely because the lossless ARDA-compatible schema adapter was used. Candidates and peptides are the canonical Parquet tables accepted by ``Evidence.from_frames`` (see ``vdjaqp/adapter.py`` and ``vdjaqp/pipeline.py``). Explicit Boolean ``background_member`` is mandatory; inserted query references must already be marked false independently of physical assay membership. No truth table or known-antibody filename reaches scoring. Source-file labels are dropped.

The model JSON contains ``conf_x``, ``conf_e``, ``substitution`` (19×19, project alphabet), ``e_db`` and optional ``n_substitutions``. Export these from an independently calibrated ``ErrorModel`` (``_conf_x``, ``_conf_e``, ``substitution``, ``e_db``, ``n_substitutions``); document its training/exclusion population privately. For held-out assessment exclude the tested antibody from calibration. This workflow never fits on the scored targets and never substitutes a technical fixture model for a scientific calibration.

Each sample produces scores, matched-peptide counts, the six existing AIRR capture tables, and a compact JSON input hash/provenance/count record. Scientific accuracy requires a separate evaluation against truth; these technical outputs make no recovery or accuracy claim.

Simulation
----------

::

   nextflow run main.nf --mode simulate-airr --repertoire /absolute/candidates.parquet --panel /absolute/panel.csv --producer 'MiXCR export' --physical_scope all --outdir results/new-simulation
   nextflow run main.nf --mode synthetic --references /absolute/raw-reference-directory --outdir results/new-synthetic

The first calls the established supplied-AIRR simulation: as-is proteins, separate proteases, MSSimulator → native search → peptide exports and existing coverage/recovery readouts. ``panel.csv`` requires exactly 12 distinct ``clone_id`` values identifying computationally inserted references (excluded from the background), not an arbitrary positive-target subset of genuinely observed background records. Physical scope (``all`` or ``panel``) never restricts the full observed search database. The synthetic technical route calls ``smoke_all.py``: GenAIRR, InSilicoSeq, fastp, **ARDA for both geometries**, then native MS/search/scoring. Generating truth stays separate from observed candidate FASTA. It is a small integration route, not the full parameter study or a new completed ARDA cohort. Install optional tools using the committed ``envs/{generation,preprocess,openms, arda}`` recipes before selecting these branches. Build the pinned ``MSSimulator-digestfix`` executable and activate the five-module ARDA overlay as documented there; defaults fail if the required native binary or map is absent. Current actual scheduled qualification is recorded in methods guide.

All modes disable task caching and reject resume and existing output roots. Use a fresh output directory; no deleted runtime is reused. ``--cpus N`` must be positive and is capped to available processors. Each process runs at most one task at a time; joint modes can run their two processes concurrently. 4 GB/15 min are scheduler declarations, **not enforced macOS memory limits**. Optional native stages also have their existing subprocess timeouts. Select a suitable executor/resource profile before large runs. The author-approved finite Aldan reproduction workload is specified in methods guide. Work directories under ``tmp/nextflow-work`` are transient. After verifying published final outputs remove that task-owned scratch; do not reuse it.

Branch dispatch
---------------

``analysis`` scores each sample-sheet row. ``real-cohort`` evaluates the supplied cohort.
``simulate-airr`` uses the supplied proteins as is. ``synthetic`` is a small tool-integration route.
``population`` runs the full synthetic population driver. ``joint`` combines population and real-cohort;
``joint-airr`` combines supplied-AIRR simulation and real-cohort. Joint modes require at least eight allocated cores,
reserve two for real analysis, and give the remainder to simulation. No automatic parameter sweep runs.

The reproduction profile declares 16 GB and 30 minutes; scientific profiler limits are 1,650 seconds,
10,000,000,000 disk bytes and 7/8 GB RSS for joint population/AIRR or 15 GB for standalone runs.
These declarations and profiler checks do not replace a scheduled allocation.

Matcher result export
---------------------

``PeptideMatcher.search`` consumes the pinned seqtree batch result's flat hit
and mismatch arrays through NumPy and Arrow. It preserves native hit order,
per-observation identifiers, substitution budgets and the existing typed placement
schema without constructing Python objects for every hit. Identical peptide
observations remain separate, including their confidence and spectrum metadata.
The search algorithm, mismatch allowance and evidence scoring are unchanged.
Optional gapped placement still lacks per-column mismatch evidence and is not
silently enabled by this export optimization.

Performance measurements for this change concern index construction, search and
placement-table materialization. They do not measure end-to-end workflow speed,
identification accuracy or sequence reconstruction.

Reproducibility qualification
-----------------------------

``workflows/qualification`` contains final-output comparison tables and source/hash
receipts. Fixed-seed qualification compares complete table inventories and values,
not merely a few headline endpoints; alternate seeds retain all existing aggregate
deltas. CPU count is part of native read-generation/acquisition provenance.
Logical decoded-read hashes distinguish compression metadata from sequence changes.
No final-output comparator is a licence to reuse generated intermediate features.

``workflows/acceptance`` retains the small software acceptance outputs. The fixture
builder ``workflows/make_fixture.py`` creates invented one-candidate/one-peptide
inputs. This acceptance is distinct from biological accuracy and native cohort
qualification. Production candidates require independently calibrated error models.

Joint execution commands
------------------------

.. code-block:: sh

   nextflow run main.nf -profile reproduction --mode real-cohort --input_root /absolute/original-input-bundle --cpus 16 --outdir results/fresh-real
   nextflow run main.nf -profile reproduction --mode population --antibodies 1200 --seed 20261110 --read_seed 2468 --cpus 16 --outdir results/fresh-population
   nextflow run main.nf -profile reproduction --mode joint --input_root /absolute/original-input-bundle --cpus 16 --outdir results/fresh-joint

Joint AIRR additionally requires ``--repertoire`` and ``--panel``. These examples
illustrate interfaces and do not authorise a new scientific run.

Declared experimental candidate populations
-------------------------------------------

``scripts/workflow_real.py --expected-candidates N`` requires exactly ``N``
candidates in every experimental context. The default remains 33,052 for the
accepted supplied-repertoire cohort. Nonpositive declarations fail before input
loading or output creation. A new raw-derived input population must be validated
and declared explicitly; this option does not authorize cohort expansion,
threshold changes or replacement of accepted conclusions. Observed variant
support, corrected clonotype membership and whole-clone abundance remain separate
quantities. Missing or ambiguous measured counts must not be replaced by zeros
or copied across compatible identities.

Native sequence-view exports
----------------------------

``scripts/export_sequence_evidence.py --input-root INPUT --output FRESH_OUTPUT
--expected-candidates N --cpus 16`` performs one fresh full-population gel-context
Evidence build and exports the existing report API's residue profiles and peptide
placements. The display selection is all twelve inserted references and four
highest full-score biological candidates, with deterministic ID ties. It is not
an identification threshold or an unbiased cohort sample. Coverage means at least
one agreeing native placement under I/L equivalence; compatible placements are
not uniquely credited evidence. The manifest declares selection and coordinate
semantics. Source data, IDs and sequences remain private.

Native parent-assignment capture
--------------------------------

``scripts/report_native_assignments.py --project-root PROJECT --output FRESH
--cpus 16`` reads the complete protected native-input and final context exports.
It collapses variant compatibility to the native heavy parent before joining
measured counts. A parent's counts are counted once per capture view, including
when several variants match. Parents without an available candidate query retain
unknown capture, with their counts reported separately. Units are native consensus
assignments and per-parent distinct-UMI assignments, not globally unique molecules,
transcript abundance or protein amount. The existing 20-nat view remains a
conditional evidence readout, without asserting unique identity.

Matched subset analyses
-----------------------

``scripts/workflow_real.py --calibration-root COMMON_INPUT`` fits each context's
error model from an explicit common input population while searching the selected
subset's candidates. The existing target-exclusion rule applies to that common
calibration too. Omitting the option preserves the original per-input route.
The receipt records common-input hashes and an exact parameter hash per context,
so a subset comparison can verify equal calibration. This does not change score
thresholds or make different search backgrounds equivalent; membership summaries
of one joint search retain their common candidate universe.

Calibration indexes only the self-reference candidates named by the retained
training observations, including database calls. The full loaded candidate
population still supplies the amino-acid composition prior. Held-out exclusions,
the substitution ceiling and repeated-offset selection remain unchanged; unrelated
background matches were discarded by the self-reference join in either route.

Separate peptide placement artifacts
------------------------------------

``scripts/workflow_real.py --mapping-output FRESH_DIRECTORY`` retains the existing
seqtree/``PeptideMatcher`` placements separately from calibrated evidence and score
aggregation. Placement rows name the peptide observation, protein-query record,
zero-based half-open interval and substitution count; a separate child table holds
substitution offsets. Query/observation indexes and hashed receipts preserve the
joins. For a scheduled single ``--case``, the supplied directory is that context's
mapping output; an all-context invocation creates one child directory per case.
Mapping is owned by these project scripts and the existing analysis pipeline,
not ARDA/MIGEC or a new standalone tool. All outputs require fresh locations.

An explicit ``--compatibility-queries PARTIAL_CATALOG`` adds a separate native
batch for observed partial intervals, after scoring memory is released. Every
query must have nonmissing false ``identity_eligible``. Its placements are stored
under ``partial_compatibility`` with a compatibility-only receipt and never enter
calibration, identity scores or acceptance. Keep the entire per-context mapping
and analysis directories within the resource monitor's output root.

Sequence-view formatting
------------------------

For the fresh raw-data DAG, use ``scripts/export_sequence_evidence.py`` with
``--retained-context CONTEXT_ROOT --display-table DISPLAY.tsv --input-root INPUT``
and the exact ``--expected-candidates`` count. The context must have complete,
hash-verified analysis and separate placement artifacts. The display table names
unique ``clone_id`` and neutral ``display_id`` values; it declares an illustrative
cohort, not an accepted-identity set or a repertoire sample. Coverage is derived
through the existing reporting API from those placements; scores are read as is.
No search, calibration or scoring is repeated. Candidate ambiguity uses the full
identity-query universe before display selection, and excludes partial-only
compatibility queries. Missing nucleotide-support metadata is labelled as a
region proxy rather than measured ancestry. The historical default Evidence-build
route remains available for archival reproduction and is not used for this DAG.

Retained calibration and alignment evidence
-------------------------------------------

Each context's provenance now retains the exact fitted error-model parameters
alongside their existing hash. Publication-mode analysis additionally writes
``alignment_evidence.parquet`` from the already computed ``Evidence.alignments``
table. This table stays in the analysis directory, separate from native placement
artifacts and aggregate scores. These exports add no search, calibration or
scoring pass and do not change scientific thresholds or scores. Historical
outputs without these exports remain historical results; their hashes alone
do not recover missing parameter values or alignment likelihoods.
