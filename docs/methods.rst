Scientific methods and interpretation
=====================================

Tool lineage and reproducibility
--------------------------------

The workflow draws on laboratory expertise and in-house MIGEC and ARDA tools to
support reproducible processing across qualified assays. MIGEC was described by
Shugay et al. (Nature Methods 2014, `doi:10.1038/nmeth.2960
<https://doi.org/10.1038/nmeth.2960>`_). ARDA is an in-house reimplementation of the
MiXCR repertoire workflow; Shugay coauthored the original MiXCR paper by Bolotin
et al. (Nature Methods 2015, `doi:10.1038/nmeth.3364
<https://doi.org/10.1038/nmeth.3364>`_). The installed ARDA release and invocation
identify the actual producer. This lineage and the pinned conventional expression
recipe describe reproducibility and versatility, without a superiority claim.


Methodological advances
=======================

The method attributes protein-derived peptide observations to candidate source clonotypes. Here “parent/source clonotype” means the repertoire sequence that could produce the protein; simulated family ancestors have a separate truth role. Probabilistic attribution is implemented through residue likelihoods, competing explanation weights and candidate-level confirmation. Preserve this scientific logic alongside the workflow and simulation machinery.

+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Advance                                          | Purpose and implementation                                                                                                                                                                                               | Evidence/status                                                                                                                                                                                                                                                                   |
+==================================================+==========================================================================================================================================================================================================================+===================================================================================================================================================================================================================================================================================+
| SHM-sensitive source attribution                 | Germline-coordinate repertoire backgrounds and error-convolved residue likelihoods weight distinguishing observed variants, including SHM-bearing positions (``background``, ``germline``, ``evidence``)                 | Implemented and integrated; evaluated as part of the joint identification method                                                                                                                                                                                                  |
+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| N-addition identity information                  | Batched junction annotation and ``pgen.insertion_pgen`` retain inferred nucleotide N-insertion lengths separately from templated D; full-junction generation probability and flank-reversion diagnostics remain distinct | Implemented standalone identity diagnostics; not silently replaced by an amino-acid core-length proxy                                                                                                                                                                             |
+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Recurrence and birthday-style collision support  | ``pgen.sharing_probability`` evaluates recurrence of a specified junction; the intended peptide-collision/birthday layer remains a separate required methodological workstream                                           | Fixed-query recurrence implemented; the selected family is all peptide–repertoire comparisons together; exclusion/null semantics remain under design review before integration. Neither empirical peptide significance nor an all-pairs collision formula is silently substituted |
+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Confidence and mass-aware error calibration      | Monotone residue-confidence calibration, a mass/composition-informed substitution matrix and separate database-call errors (``errormodel``)                                                                              | Implemented; target-excluded evaluation retains eight individual runs without fitting their target readings                                                                                                                                                                       |
+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Shared-observation competition                   | Candidate/offset explanations compete with each other and a null explanation (``aggregate``)                                                                                                                             | Implemented and represented in existing sequential method comparisons                                                                                                                                                                                                             |
+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Bounded residue confirmation and contradiction   | Agreeing coverage uses a soft union and per-residue information cap; contradictory observations contribute bounded negative evidence (``aggregate``)                                                                     | Implemented; prevents repeated overlapping peptides from being counted independently in full                                                                                                                                                                                      |
+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Isotype prior and junction-information reporting | Group-normalized IgG prior, separate confirmed/junction nats and a reporting-only framework flag (``aggregate``)                                                                                                         | Implemented; keep its contribution and preparation-specific trade-offs visible                                                                                                                                                                                                    |
+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Explicit observable identity                     | I/L-equivalent full-variable identity classes, partial supported contigs and separate physical/observable/recovered denominators (``contigs``, ``simulation_metrics``)                                                   | Implemented; avoids inferring complete identity from partial compatibility                                                                                                                                                                                                        |
+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Native known-truth and versioned RNA workflow    | Published repertoire/read/MS tools, related SHM families, peptide-contributor truth and observed boundary handling                                                                                                       | Original twelve coding-only and six leader-augmented RNA builds complete, with physical-payload invariants; archival TRUST4-consensus/ARDA-annotation producer, distinct from current ARDA-only qualification                                                                     |
+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Experiment-conditioned reconstruction            | Supplied protein sequences retained as is; explicit physical/search membership; yield and per-reference coverage drive bounded pyABC inference                                                                           | Forward model qualified; inference status belongs to the active execution record. Recovery and newly reported design fractions are excluded from fitting                                                                                                                          |
+--------------------------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Definitions and preservation rules
----------------------------------

SHM attribution uses observed sequence evidence; the final aggregate score is not claimed to be a calibrated whole-clonotype posterior. The current N-length helper counts nucleotide slices of an inferred Viterbi/scenario reconstruction from an amino-acid junction. Its 4^(-n) heuristic is an ordering proxy, not a calibrated peptide probability integrating codon degeneracy. Supplied junction-core label 6 can include templated D and must not be called pure N insertion. A specified-query recurrence probability and an all-pairs birthday probability count different comparison events. The intended collision layer must name that event, search population and assumptions before a “definitive match” interpretation.

Optional and standalone components are distinguished from integrated scoring. Raw experimental data and source mappings are not distributed with the software.

Contextual N-support design update
----------------------------------

The selected readout requires a minimum run of two or three consecutive agreeing N-derived amino acids within an existing contextual whole-peptide placement. Reuse ``PeptideMatcher.search`` and pinned seqtree 1.1.0 batch results, including every relevant offset and mismatch positions; no short-motif-only index is needed. Three mapped-peptide roles are distinct: previously accepted evidence, germline-explainable evidence and exclusion from null estimation. Treat those roles separately, retain existing main results, and define inferred nucleotide/codon footprints without relabelling templated D as N. Numerical collision significance is not supplied by this contextual support predicate.

Preserve source background membership across identifier conversion
==================================================================

Candidate population role is independent of candidate names. Own-residue subtraction applies only to valid residue slots that actually contributed to the background counts.

Root cause and interface
------------------------

``Background`` now accepts an optional candidate column **``background_member``**, whose type must be Boolean and whose values must be nonmissing. When present, it supplies the row membership directly. When absent, the existing ``exclude_prefix`` behavior is unchanged. Own-residue subtraction still occurs only for valid residue slots that actually contributed. No frequency, smoothing, coordinate, error model, prior, aggregation or threshold formula changes.

``vdjaqp.airr_ms_inputs.attach_background_membership(candidates, inserted_reference_ids=...)`` attaches the role from the explicit source reference inventory before identifier replacement. References must be unique and present; existing explicit membership is preserved, with an error if an inserted reference is marked as a member. The helper does not inspect which proteins will be physically simulated. The twelve inserted-reference source IDs in the AIRR callers define the excluded role even when the physical input contains the entire repertoire. Physical membership and background membership remain separate columns.

The candidate snapshot retains the added Boolean column; the private identifier mapping carries it through the bijection, and the input manifest records member/excluded row counts. Existing source sequence, nucleotide and annotation fields, original source IDs, row order and candidate membership are not replaced by this operation. Search/protein FASTA sequences and intensity allocation use the same existing fields and rules. Source target-file hashes remain in the joint provenance; the added manifest count is not a substitute for that source record.

Caller and adapter trace
------------------------

+--------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Boundary                                                                       | Behavior                                                                                                                                                                                                                                                                                                                  |
+================================================================================+===========================================================================================================================================================================================================================================================================================================================+
| ``vdjaqp.pipeline.Evidence.from_frames``                                       | Passes the candidate frame to the shared ``Background``; no special scoring bypass or per-caller probability formula.                                                                                                                                                                                                     |
+--------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``vdjaqp.evidence.AlignmentScorer``                                            | Uses an explicitly supplied ``Background``, or constructs that same shared implementation.                                                                                                                                                                                                                                |
+--------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``vdjaqp.abundance_forward.forward_draw``                                      | Attaches source reference roles before writing inputs or neutralizing IDs. All three acquisition readouts use the role-preserving candidate table in subsequent scoring.                                                                                                                                                  |
+--------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``scripts/simulate_airr_ms.py``                                                | Attaches source reference roles before physical-scope selection and identifier conversion. Both panel and whole-repertoire physical branches retain the same excluded source-reference roles.                                                                                                                             |
+--------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| ``prepare_airr_ms_inputs`` / AIRR-compatible schema adapter                    | Preserve extra candidate columns; serialize the explicit role and copy it into the private ID mapping when present. Physical selection never determines this role.                                                                                                                                                        |
+--------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Generic observed GenAIRR scoring via ``smoke_score.py`` / ``smoke_metrics.py`` | No role is inferred from known truth, simulated abundance or the physical assay panel. Existing genuine generated/observed candidates retain their existing inclusion behavior.                                                                                                                                           |
+--------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Legacy experimental ``Evidence.build``                                         | No new column is required; original source-ID prefix handling remains unchanged.                                                                                                                                                                                                                                          |
+--------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Standalone ``reconstruction.append_known_references``                          | Genuine input candidates receive true when the optional role is absent; existing explicit roles are validated and preserved. Appended supplied references receive false independently of their IDs. The original ``candidate_construction`` provenance remains, and all pre-existing input payload columns remain intact. |
+--------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Selectable aggregation of competing spectrum explanations
=========================================================

``Evidence.score(aggregation="exclusive_spectrum")`` selects the alternative described below. ``aggregation="alignment_union"`` remains the default. This is an explicit model alternative, not a claimed improvement in precision, recall or ranking. Promotion requires matched scientific evaluation; constructed probability examples and software timings cannot establish that improvement.

Scope and equations
-------------------

For an alignment a, let s(a) be its spectrum ID, c(a) its candidate, and γₐ its unchanged spectrum-level responsibility. The existing responsibility calculation includes a null component and competes all returned readings, candidates and placements within s (``vdjaqp/aggregate.py:165``). Let Iₐᵣ equal one when alignment a covers candidate residue r with an agreeing canonical encoded amino acid, and zero otherwise. Candidate residue slots are globally unique across the concatenated candidate sequences. Existing encoding includes I/L equivalence; unknown/negative codes never agree.

The historical heuristic treats each agreeing alignment as a separate union event:

::

   C_old(c,r) = 1 − ∏_{a:c(a)=c} (1 − γₐ Iₐᵣ).                 (1)

The alternative first sums mutually exclusive explanations within the same spectrum, then unions across distinct spectrum IDs:

::

   q(s,c,r) = Σ_{a:s(a)=s,c(a)=c} γₐ Iₐᵣ,                   (2)
   C_new(c,r) = 1 − ∏_s (1 − q(s,c,r)).                      (3)

All γ, q and C are dimensionless model quantities. The existing null-inclusive competition gives Σ\_{a:s(a)=s}γₐ ≤1, so q is in [0,1] up to floating-point error. The implementation retains the existing ``max(1−q, 10⁻¹²)`` log-miss floor, including its endpoint approximation. A spectrum with no agreeing alternative contributes zero q; an absent spectrum/residue pair contributes a multiplicative one. Empty scored input gives zero confirmation.

Two same-spectrum alternatives with γ=.4 each give q=.8; two distinct spectra with γ=.4 each give C=.64. These illustrate the change in assumed event structure, not biological accuracy or posterior calibration. Product combination across distinct spectrum IDs is a model independence assumption; distinct keys alone do not prove experimental independence. Shared peptides still do not establish a unique candidate identity.

Contradiction remains the existing capped linear sum of γₐ min(LLRₐ,0) per candidate. Per-residue information weights/caps, junction flags, isotype prior, final ranking/ties and all recovery thresholds are unchanged (``aggregate.py:268`` onward). No SHM or junction bonus is added. Scores and coverage may change because the responsibilities assigned to alternatives are combined differently; the raw responsibilities themselves do not change.

No-competition semantics and callers
------------------------------------

With ``competition=False``, responsibilities are replaced by ones. Those are deterministic coverage indicators, not a categorical distribution of competing explanations. Both aggregation selections therefore explicitly use the historical alignment-union path for this ablation. This preserves its behavior instead of summing ones as though they were probabilities.

``Evidence.score(**kw)`` already forwards these arguments to ``ClonotypeScorer.score``; no pipeline interface change is required (``vdjaqp/pipeline.py:102``). Existing benchmark, report, simulator and abundance-forward callers omit ``aggregation``, so retain the historical default. A comparison must explicitly request the alternative and record the selected value. Explicit native ``spectrum_id`` values remain preferred and validated; when absent, the existing file/mass/charge/RT grouping heuristic is unchanged (``aggregate.py:64``). This change does not redefine observation keys or assert that legacy heuristic groups are known physical spectra.

Contextual inferred-N support reducer
=====================================

``vdjaqp.junction_support.contextual_n_support`` consumes existing whole-peptide placements; it performs no motif search, native annotation, Pgen calculation, fitting, scoring or empirical-null generation. It does not modify accepted scores. The minimum consecutive AA run is an explicit parameter, default three and configurable to any positive integer (the initial comparison includes two and three).

Inputs and coordinates
----------------------

All arguments are Polars tables. Unique catalogues load sequence strings once; placement rows do not duplicate peptide strings.

+-----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Table                             | Required columns and meaning                                                                                                                                                                                                                                       |
+===================================+====================================================================================================================================================================================================================================================================+
| candidates                        | ``clone_id``, canonical protein ``seq``, ``cdr3_start`` (zero-based AA offset in that exact candidate), ``cdr3_seq`` (anchor-inclusive observed junction). Absent junction fields can remain null and produce unavailable geometry, including germline competitors |
+-----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| peptides                          | Unique ``peptide_id``, canonical AA ``seq``                                                                                                                                                                                                                        |
+-----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| placements                        | ``observation_id``, ``peptide_id``, ``clone_id``, zero-based half-open AA ``start,end``, Boolean ``is_germline``. Length must equal the entire peptide; no isolated-run query is accepted                                                                          |
+-----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| annotations                       | ``clone_id``, inferred ``cdr3_nt``, native ``v_end_nt,j_start_nt,d_start_nt,d_end_nt,np1,np2``                                                                                                                                                                     |
+-----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| accepted_scores                   | Unique ``clone_id``, finite ``confirmed_nats`` from the existing accepted scoring output                                                                                                                                                                           |
+-----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| accepted_support                  | ``observation_id,clone_id``: the caller’s existing **credited-support relation**, not arbitrary compatible alignments. The reducer cannot infer credited support from a placement alone                                                                            |
+-----------------------------------+--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Observed candidate offsets come from existing candidate records, whose anchored ``cdr3_seq`` must exactly match ``seq[cdr3_start:cdr3_start+len(cdr3_seq)]``. They are never supplied from simulation truth or guessed from germline. The annotation must correspond to that exact AA junction; inferred NT translation is validated. Inferred nucleotide reconstruction remains labelled inferred, not measured N bases or a complete posterior over boundary scenarios.

The adapter follows the inspected installed ``vdjtools.model.junction.annotate_junctions`` return construction exactly: first N tract is ``cdr3_nt[v_end_nt:d_start_nt-1]``, second is ``cdr3_nt[d_end_nt:j_start_nt]``. D coordinates are one-based closed; the resulting N slices are zero-based half-open. Both slices must equal the native ``np1,np2`` strings (null permits a genuinely empty slice only when complete coordinates validate). Missing NT or any required coordinate yields ``missing_geometry``; it is never interpreted as a zero-length N assignment. Inconsistent nonmissing geometry raises before support is produced.

No native annotator is called by this module. Callers can obtain annotations through the existing **single batched** ``annotate_junctions`` API, after its provenance/sequence correspondence is qualified. This implementation does not read annotations from a computed-result cache.

Eligibility and retained roles
------------------------------

Each AA is eligible only if its whole codon lies inside **one** inferred N tract. Mixed N/templated codons and D-overlapping codons do not contribute. Candidate and peptide positions must agree under the existing I/L-collapse alphabet. Mismatches, non-N positions and D/mixed positions break consecutive runs; discontiguous agreement cannot be added to meet the minimum.

Every input placement remains represented in input order with an explicit placement index, identifiers, interval, germline label, geometry status/provenance, longest eligible run and counts of agreeing N, mixed and D-overlapping residues. Germline explanations remain competitors/labels and never become rescue discoveries. An observation credited to any existing candidate with confirmed information ≥20 nats is masked globally from this rescue view. Original placements, accepted scores and their membership are unchanged; the module does not remove ordinary scorer explanations or choose a new threshold.

``rescue_eligible`` is a Boolean eligibility view, not an attribution probability, significance value, new accepted identity or altered score. Future family-wide collision statistics and empirical-null modelling remain unimplemented. The three roles—explained-observation mask, contextual support predicate and retained germline explanation—are kept distinct.

Exact reference-backed IGH V annotation
=======================================

``generate_population.py`` rebuilds an allele-level map once per invocation from installed raw GenAIRR references. It verifies both SHA256 pins in ``envs/generation/identity.json`` before creating the output directory or loading the trusted legacy pickle. No generated sequences, native configuration, draws, gene weights, or native ``airr_truth.v_call`` values change.

The source is ``builtin_dataconfigs/HUMAN_IGH_OGRDB.pkl`` (198 V alleles); the target is ``published_references/human_igh_ogrdb_10_observed_v1.zip``, member ``reference.genairr-ref.json`` (190 V alleles). ``source.lock.json`` is retained in truth provenance. One in-memory SQLite left join compares uppercase ungapped nucleotide strings exactly. No alignment, inferred group rename, cross-run cache, or reference switch is involved. Legacy alleles from a shared native gene group can map to different canonical genes.

Each chain gains three scalar fields:

+-----------------------------------+-------------------------------------------------------------------------------------------------+
| Field                             | Meaning                                                                                         |
+===================================+=================================================================================================+
| ``canonical_v_gene``              | Canonical gene string only when exactly one distinct gene matches; JSON null otherwise          |
+-----------------------------------+-------------------------------------------------------------------------------------------------+
| ``canonical_v_match_status``      | ``unique_gene``, ``ambiguous_genes``, ``no_exact_match``, or ``not_evaluated_locus``            |
+-----------------------------------+-------------------------------------------------------------------------------------------------+
| ``canonical_v_gene_count``        | Count of distinct exact matching genes: 1, greater than 1, or 0; JSON null for unevaluated loci |
+-----------------------------------+-------------------------------------------------------------------------------------------------+

IGH is evaluated. IGK and IGL have status ``not_evaluated_locus``, not a claimed absence of mapping. Multiple allele matches to one gene remain ``unique_gene``. An unknown IGH native allele fails explicitly. In the pinned metadata, among 198 legacy V alleles, 163 map to one gene, seven map to multiple genes and 28 have no exact match. These are reference metadata counts, not population yields.

Top-level ``truth.json.reference_gene_mapping`` records the method, locus/segment, both reference hashes, source lock, reference allele counts and all 198 legacy alleles keyed by native allele name. Each entry preserves the legacy gene, sequence SHA256, complete canonical allele and distinct gene lists, explicit allele/gene pairs and scalar annotation. No nested mapping is added to a chain. ``reference_gene_mapping.tsv`` is the flat lookup deliverable: one exact pair per row, including an explicit ``null`` canonical allele/gene row for every unmatched legacy allele. Its columns are ``legacy_allele``, ``legacy_gene``, ``sequence_sha256``, ``canonical_allele``, ``canonical_gene``, ``canonical_v_match_status`` and ``canonical_v_gene_count``. Native calls remain the join key for chain truth.

.. code:: sh

   .tool-envs/generation/bin/python -m unittest discover -s tests -p test_genairr_gene_map.py -v

The pure fixture join test also runs without GenAIRR in the baseline interpreter; the installed reference check is skipped there when GenAIRR is absent.

The evaluation metadata and gene-coverage exports carry these three scalar fields while retaining native ``v_call``/``gene_call`` keys. Older fixtures without a mapping are ``not_evaluated``. D/J/C summary rows use ``not_applicable`` for this V-only annotation. Summaries keep their original native bins, chain counts and residue denominators; ambiguous pairs are not exploded into duplicated observations. The flat lookup uses the literal ``null`` for unmatched canonical names, accompanied by ``no_exact_match`` and count zero. A uniquely named reference match is not proof of genomic locus identity or productive generation.

Conditional paired population generator
=======================================

``scripts/generate_population.py`` reuses the validated ``smoke_generate.load_constants`` and ``reference_hashes`` APIs without modifying the historical fixture. Run with the pinned generation interpreter:

.. code:: sh

   .tool-envs/generation/bin/python scripts/generate_population.py \
     --output /absolute/fresh/absent-directory --reference-dir /absolute/raw-six-records \
     --families 40 --pairs-per-chain 100 --mutation-rate .01 \
     --protein-concentration equal --seed 20261110
   .tool-envs/generation/bin/python tests/test_population_generator.py

The callable is ``generate(output, reference_dir, families, pairs_per_chain=100, mutation_rate=.01, protein_concentration='equal', seed=20261110, heavy_constant_gene='IGHG1', additional_reference_dir=None)``. Output must not already exist. Counts must be positive integers, families even; rates finite in [0,1]. Raw references must match all six exact hashes from the existing generation manifests before any output directory is created. For the default, only these six raw records are copied. Every invocation builds all native outcomes afresh.

Each family comprises three paired antibody identities, root → child → grandchild, with heavy and light sharing ancestry. Half the families use kappa and half lambda, an assumed balance. By default all physical heavy constants are selected IgG1; lambda uses the existing native ambiguous J2/3 allele pool restricted to IGLC2. This is a conditional physical IgG1 population, not unrestricted genomic J–C sampling or the full matched RNA population. Forty families imply 120 antibody identities, 240 chains and 24,000 pairs per geometry at 100 pairs per chain.

An explicit ``--heavy-constant-gene IGHM|IGHG3|IGHA1|IGHE`` selects one fixed existing secreted reference for every heavy chain. Supply ``--additional-reference-dir /absolute/pinned-public-raw-inputs`` with all raw files listed in ``envs/generation/constant_candidates.json``. The existing ``validate_constant_candidates.load_validated_constants`` validates every source hash, annotation, translation and join before output creation or native calls. An unknown gene, missing explicit directory, or an additional directory with default IGHG1 fails. The default still uses the unchanged six-record loader, seeds and draw policy.

Additional raw files are copied under ``references/additional_constants``, alongside the exact source manifest; truth records their hashes, source URLs, chosen gene, and annotated CDS/protein identities. IGHM/IGHG3 overlap one identical native terminal G (the appended suffix excludes the duplicate G); IGHA1/IGHE complete the split initial codon. Whole-draw compatibility checking applies without base repair. IGHG3 retains its source’s natural three-exon hinge variant; IGHE retains its source’s curated W43C variant, without substituting the displayed UniProt sequence. These choices specify fixed reference sequences, not class frequencies, alternative isoform sampling, mixed-class populations or polymerization assumptions.

Roots use one native compiled ``run(n, seed, strict=True)`` batch per locus. Native batch seeds are ``seed + locus_index*10^14 + batch_index``. Retries and child continuations use injective mixed-radix stage/locus/family/node/attempt coordinates; all coordinates and accepted seeds are recorded. The declared family bound keeps those ranges disjoint. A child always continues its accepted native parent through ``simulator.run_from``; heavy/light seeds are independent. No homemade recombination, mutation, affinity, errors or digestion are used. Current pinned native batch uses a serial map, not Rayon; no process wrapper is added. Profile end-to-end before introducing a parallel layer.

Full native variable nucleotide sequence plus annotated extracted C CDS must translate exactly to native variable amino acids plus the fixed reference constant and terminal stop. Native productivity alone does not guarantee this split-codon boundary. An incompatible **whole draw** is rejected; another native draw uses the same parent and next predeclared attempt seed. Roots are conditioned likewise. No base repair, reverse translation, J masking or genotype omission occurs. Each node has at most 100 attempts, then generation fails explicitly and preserves ``failed_attempts.json``. Unrelated engine errors propagate. This acceptance conditioning changes the unconditional native distribution; rejected terminal mutations are not asserted biologically impossible. All native attempts, reasons, seeds and raw mutation counters remain in truth. Rejection counts remain explicit in each run’s truth records.

``airr_truth`` preserves native counters verbatim. Separate typed branch/root distance objects contain actual equal-length nucleotide Hamming substitutions, denominator nt and fraction, including the terminal native nucleotide. Native observable mutation counters are never relabelled as root distance. The model permits substitutions only; length changes fail explicitly. AIRR junction includes conserved anchors; no alternate CDR3 extraction is invented.

RNA pairs are externally prescribed equally for every chain in both protein conditions. Equal protein concentration is 1,000 arbitrary molar intensity units per chain. ``head80`` gives the first ceil(0.1 × families) families in ascending family ID 80% of the same total pool and the remainder 20%, equally within each family. ``rna_pairs * secretion_factor == protein_intensity``; these factors are imposed, not measured secretion or affinity. Paired H/L weights are equal; a common 2H2L factor is uncalibrated. The concentration perturbation is a mechanistic assumption, not a fit to experimental abundance.

``constant_reference_proteins.fasta`` stages the three selected public constant proteins as declared search model inputs, before native generation. Its sequences come only from validated reference records, never generated variable-chain truth; headers are ``constant_{protein_accession}``. Truth names this artifact. The ``smoke_metrics.py --families ...`` driver forwards both constant options and uses this staged reference FASTA for both search databases; a missing file fails instead of reverting to IgG1. Its legacy fixture route retains its original three constant-reference records. No generated variable truth enters either database.

Outputs retain ``truth.json`` chain keys, ``proteins.fasta`` intensity headers, ``amplicon_templates.fasta``, ``rnaseq_templates.fasta`` and ``read_counts.tsv`` (twice RNA pairs, the existing read generator convention). Truth adds parameters, constant/source hashes, pinned source identity, runtime versions, seed policy, conditioning records and actual ancestor distances. Coding RNA includes terminal stop; no UTR, signal peptide, splicing efficiency or assembly model is added here.

Explicit naive/memory mixture
-----------------------------

Omit ``--families`` and supply ``--total-antibodies 1200 --memory-fraction .25`` (or ``.5``/``.75``) plus the pinned ``--additional-reference-dir``. The callable accepts these same named arguments; mixing them with ``families`` fails before writing. Mixtures use fixed IGHM for naive pairs and IGHG1 for memory pairs. Native recombination produces independent naive singleton pairs without mutation. Each memory family has a latent unmutated ancestor followed by three observed native S5F continuations in a chain. Latent ancestors remain in truth only and never enter proteins, templates, read counts, or the observed antibody denominator. Neither component promises sequence uniqueness after independent rearrangement draws.

+-----------------------+----------------------+-----------------------+------------------+
| Memory count fraction | Naive observed pairs | Memory observed pairs | Memory families  |
+=======================+======================+=======================+==================+
| 0.25                  | 900                  | 300                   | 100              |
+-----------------------+----------------------+-----------------------+------------------+
| 0.50                  | 600                  | 600                   | 200              |
+-----------------------+----------------------+-----------------------+------------------+
| 0.75                  | 300                  | 900                   | 300              |
+-----------------------+----------------------+-----------------------+------------------+

These counts always total 1,200 observed antibody pairs and 2,400 observed chains. Kappa/lambda are balanced separately within each component. Native Poisson SHM can produce zero realized substitutions; component labels do not replace measured branch/ancestor distances. Memory node zero is already an observed continuation, not the unmutated latent ancestor. Chain truth adds ``component``, ``ancestor_id`` and ``is_head``; ``parent_id`` can refer to a latent ancestor. ``latent_ancestors`` and ``component_counts`` record the separate denominators. The root-distance origin for memory is its latent ancestor; naive singleton distances are zero.

``head80`` selects exactly ten percent of independent naive singletons and ten percent of whole memory families, separately within each light locus. Heads carry 80% of protein mass inside each component/light stratum, retaining each stratum’s total and the requested memory fraction. This requires each component family count to be divisible by 20; incompatible sizes fail rather than round. Equal-concentration fixtures only require even naive counts and memory counts divisible by six. RNA pairs remain equal per observed chain in both abundance settings. The declared constant model-input FASTA contains both heavy references and both light references. Root seeds retain locus-separated native batches; latent-root retries use stage one, and three observed memory continuations use stage two with nodes zero through two.

Observation-derived transcript boundaries
=========================================

``vdjaqp.transcript_boundary.project_transcript_boundaries(queries, contexts, alignments)`` projects verified public mature-start markers through supplied observed alignments using ``arda._markup.transfer_regions``. It performs no search, genotype selection, threshold filtering, sequence completion or adapter cropping.

Queries provide ``sequence_id``, uppercase ungapped ``sequence``, and ``retained_primary_exact_match_mask``. Sequence and Boolean mask must already share the coding orientation; missing masks leave boundaries unresolved. Contexts provide ``target``, ``sequence``, and integer ``mature_start`` in that exact public context. Upstream preparation owns splicing, physical mature-start interpretation, source offsets and any difference from ARDA’s frame-trimmed V scaffold.

Alignment rows provide ``sequence_id``, ``target``, integer ``qstart``, ``qend``, ``tstart``, ``tend``, ``qlen``, ``tlen``, and ``qaln``, ``taln``. Coordinates are nucleotide, 1-based closed, on forward coding-oriented sequences. Ungapped alignment strings must equal their source slices. All returned retained alternatives must be supplied, including nonbracketing alternatives; this component cannot verify a caller’s search limits or omitted hits.

Both target bases flanking the boundary must project to adjacent query bases, and both query bases must be canonical ACGT with retained-primary exact-match support. Query N is allowed elsewhere; public contexts must be unambiguous ACGT. Every retained alternative must establish the same boundary. The flat result retains no-hit queries and reports nullable ``mature_start``, ``status``, alignment count, distinct target candidate count, bracketing/nonbracketing alignment counts and supported alignment count. Status is ``resolved``, ``no_alignment``, ``unsupported_boundary`` or ``conflicting_boundary``. Conflicting supported boundaries take precedence over unsupported alternatives. Aliases agreeing on the boundary are retained and counted.

A resolved nucleotide boundary does not certify codon phase, downstream residue support or biological cleavage. The ordinary contig adapter is unchanged.

The optional ``search_transcript_boundaries`` wrapper makes one native batch against all 547 model contexts: donor prefix plus each native allele’s own core. It takes only observed raw contigs, observed annotation/orientation, read-support masks and the global model catalogue. It accepts no generated-chain truth or genotype selection. It uses ARDA’s existing highest whole-target bit-score ranking (``arda.annotate.mapper._best_hits``), preserving every exact-maximum tie rather than ARDA’s single deterministic winner. The pure projection helper still requires consensus across every supplied retained alternative. Weaker returned hits do not veto a better-supported boundary. Returned, retained (``alignment_count``) and lower-score alignment counts remain explicit alongside ``best_bits``; requested raw alignment diagnostics are gzip-compressed. The existing max-seqs 300 filter bounds returned alternatives; retained consensus is not exhaustive across every possible target/alignment. No scoring threshold or sweep was added.

``adapt_contigs(..., boundaries=frame)`` optionally excludes upstream leader labels and extends a supported FR1 prefix only to a resolved observed boundary. Its existing whole-codon translation and exact-primary masks still govern every residue. Missing/unsupported/conflicting boundaries preserve the ordinary path; they do not prove leader absence or complete mature-variable observability. Omitting ``boundaries`` preserves the existing adapter behavior.

Natural human constant-reference support
========================================

Existing library route
----------------------

Installed ARDA ``refbuild/constant.py`` documents ``database/c_genes/<organism>.fasta`` as **CH1 exons**, used for J–C read scaffolds (150-nt stubs). Its ``germline.segment_germlines`` supports V and J, not full spliced constant transcripts. These APIs therefore cannot supply the required full physical constants. No new sequence assembly engine is needed: Biopython feature extraction and translation suffice for integrity checks.

Verified natural candidates
---------------------------

All four reviewed UniProt entries display **membrane isoform 2**. Their explicit ``-1`` isoforms are secreted. Never use bare-accession FASTA as if it were the secreted constant. Protein/CDS checks below concern constant fragments, without variable regions, leader peptides or UTRs.

+-----------------+--------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Gene            | Secreted protein                                                         | Natural nucleotide record / annotated protein                                                                                                                 | Verified integrity and remaining use restriction                                                                                                                                                                                                                                                                                                                                                                                                                        |
+=================+==========================================================================+===============================================================================================================================================================+=========================================================================================================================================================================================================================================================================================================================================================================================================================================================================+
| IGHM            | `P01871-1 <https://rest.uniprot.org/uniprotkb/P01871-1.fasta>`__, 453 aa | `ENST00000390559.6 / ENSP00000375001.2 <https://rest.ensembl.org/sequence/id/ENST00000390559?type=cds&content-type=text%2Fx-fasta>`__, 1362 nt including stop | Exact curated secreted UniProt protein. Four annotated GRCh38 CDS features independently extracted from genomic sequence by Biopython CompoundLocation reproduce the entire downloaded CDS exactly. Native identical-G overlap passes complete variable-plus-C translation. Selected reference.                                                                                                                                                                         |
+-----------------+--------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| IGHG3           | Natural source protein CAC10261.1, 362 aa                                | `AJ390279.1 / CAC10261.1 <https://www.ncbi.nlm.nih.gov/nuccore/AJ390279.1>`__, six-part annotated CDS, 1089 nt including stop                                 | Complete natural **three-exon hinge** constant, documented in the source study; not the four-exon hinge UniProt displayed variant. Exact annotated and full native variable-plus-C translation pass with one-base G overlap. Selected candidate under this exact accession; no IMGT allele inferred.                                                                                                                                                                    |
+-----------------+--------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| IGHA1           | `P01876-1 <https://rest.uniprot.org/uniprotkb/P01876-1.fasta>`__, 353 aa | `J00220.1 / AAC82528.1 <https://www.ncbi.nlm.nih.gov/nuccore/J00220.1>`__, three-part CDS, 1061 nt including stop                                             | ``codon_start=3``; extracted CDS translates to UniProt residues 2–353 plus stop. Native terminal J G + complete extracted CDS translates **exactly** to all 353 UniProt residues plus stop. This is the verified secreted split-codon candidate.                                                                                                                                                                                                                        |
+-----------------+--------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| IGHE            | `P01854-1 <https://rest.uniprot.org/uniprotkb/P01854-1.fasta>`__, 428 aa | `J00222.1 / AAB59395.1 <https://www.ncbi.nlm.nih.gov/nuccore/J00222.1>`__, four-part CDS, 1286 nt including stop                                              | ``codon_start=3``; exact annotated translation. Native terminal J G completes first A, but natural record has C at constant residue 43 whereas displayed UniProt secreted sequence has W. UniProt VAR_044229 explicitly identifies W43C as \**IGHE*01*\*, supported by the cited genomic study. Select this exact natural source allele: 428-aa physical constant, with full variable-plus-C translation verified. Do not label its protein sequence as exact P01854-1. |
+-----------------+--------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

UniProt metadata pins: IGHM sequence version 5 (entry version 225); IGHG3 sequence version 3 (entry 207); IGHA1 sequence version 3 (entry 226); IGHE sequence version 2 (entry 205). The downloaded records and manifest checksums are authoritative if database metadata changes. Identical protein sequence does not by itself prove nucleotide allele identity; no new IMGT allele assignment is asserted.

Reference checks verify one native retained IGH terminal G convention, not all J alleles or genomic combinations. For ``codon_start=1``, the reference includes the complete first C codon: assert its first G equals the native terminal J G and join ``coding_variable_nt + constant_CDS``, equivalently ``variable_nt + constant_CDS[1:]``. This is an identical-base overlap, not an allele patch or duplication. For ``codon_start=3``, retain the complete extracted CDS and join ``variable_nt + constant_CDS``. Both routes are checked against the native variable amino acids plus the **actual source-specific physical constant** plus its annotated terminal stop. Full cohort generation must validate every chosen native J tail and complete coding translation. Native C-call labels alone do not resolve membrane versus secreted RNA isoforms. RNA background transcripts and secreted physical protein pools need separately declared isoform choices; physical-pool enrichment does not imply the same RNA isotype distribution.

Exact retrieval and integrity check
-----------------------------------

Raw official records and the reference-only assertion script are in ignored ``tmp/constant-reference-audit/``. Commands run from the repository root in the pinned generation environment:

The script uses published Biopython ``SeqIO``, ``SeqFeature.extract`` and annotated ``codon_start``, asserts every selected annotated translation, tests exact J-G/CDS/secreted-protein agreement and records raw-file/CDS SHA256 checksums. No reverse translation, allele repair, stop addition or isoform mixing occurred. ENA timeouts were transport failures; successful NCBI records are the official nucleotide source used here.

Unresolved choices and limits
-----------------------------

IGHM source closure
-------------------

`UniProt P01871 <https://rest.uniprot.org/uniprotkb/P01871.txt>`__ explicitly links ENST00000390559.6 and ENSP00000375001.2 to secreted isoform P01871-1. The `Ensembl transcript annotation <https://rest.ensembl.org/lookup/id/ENST00000390559?content-type=application%2Fjson>`__ identifies IGHM-201, HAVANA source, GRCh38, reverse strand. Its four `CDS features <https://rest.ensembl.org/overlap/id/ENST00000390559?feature=cds&content-type=application%2Fjson>`__, in transcript order, are chromosome 14 positions 105855906–105856218, 105855480–105855815, 105854917–105855234 and 105854343–105854737 (1-based inclusive, all reverse strand). Their lengths are 313, 336, 318 and 395 nt, respectively, summing to 1362 nt. Biopython converts these to zero-based half-open reverse-strand CompoundLocation parts; extraction exactly equals the downloaded versioned CDS. Translation is exactly all 453 P01871-1 residues plus stop. No inferred bases, edits or synthetic stop were introduced. These features specify this secreted transcript; membrane isoforms are separate and are not substituted.

All four requested heavy gene classes now have verified natural secreted references and full coding-boundary translation assertions. This closes reference integrity for a declared secreted candidate-RNA model; it does not establish empirical isotype frequencies or validate all possible native J/allele combinations. IGHG3 and IGHE retain their explicit source-specific natural variants. Native repertoire generation and cohort-level acceptance checks remain future work.

Native terminal nucleotide during SHM continuation
--------------------------------------------------

The native terminal G is coding sequence at a genuine split J–C codon, not removable padding. Official `AH005277.2 <https://www.ncbi.nlm.nih.gov/nuccore/AH005277.2>`__ annotates its J4 CDS at 1-based ``<1914..>1959`` with 14 complete translated residues and one remaining nucleotide G, immediately followed by the annotated J4–C intron at ``1960..>3182``. The local genomic boundary is ``...CCGTCTCCTCAG|GTGAGTCCTCA...``: the final G is exonic, immediately before the donor’s intronic GT. Together with the official C CDS features above, this supports the reference split-codon convention. It does **not** establish that every substitution of that exonic G is biologically impossible; these sequence-generation APIs do not model splice efficiency or transcript survival.

For the declared unchanged natural C-protein model, preserving G is necessary: the first constant codon is GGG for selected IGHM and GCT/GCA/GCC for selected IGHG3/IGHA1/IGHE. Replacing that shared first base with A, C or T changes the first constant amino acid, while leaving the variable coding frame and its complete translated codons unchanged. Those alternative first codons are not stops. Consequently a variable-region productivity check cannot establish exact constant-protein identity or splicing competence. If a changed G is accepted as an actual biological mutation, its altered first constant residue must be represented consistently in full RNA and protein truth; it cannot be silently replaced by the reference G or presented as the unchanged reference constant. That expanded model would require a separately justified splice/translation policy.

Installed GenAIRR 2.4.0, pinned source commit ``ae6e7d8e56bf593a87718bf493606dcb703afe87``, was inspected without generation. Its published ``Experiment.mutate`` signature has ``segment_rates`` and ``v_subregion_rates``, but no nucleotide-position mask. ``_validate_segment_rates`` accepts the V/D/J/NP buckets with nonnegative multipliers, and ``_MutateStep`` documents their application to those biological segments. Setting ``segment_rates={"J": 0.0}`` is the existing native route to exclude **all** J sites from new SHM; it is a broader modeling assumption than protecting one terminal site and must not be silently substituted. It also does not repair a parent state whose terminal base was already altered. V-subregion rates do not address J sites.

``productive_only()`` documents V/J anchor residue classes, J-tail **frame** and absence of stops in complete observed coding codons. It explicitly does not predict expression or receptor function; it is not an exact terminal-G constraint. Its Python builder exposes only the productive constraint bundle; the installed native module exposes ``ContractSet`` and the ``productive`` factory, without a documented arbitrary-position preservation factory. Source: the `published GenAIRR repository at the pinned commit <https://github.com/MuteJester/GenAIRR/tree/ae6e7d8e56bf593a87718bf493606dcb703afe87>`__; installed ``GenAIRR/_experiment/mutation.py``, ``_experiment/constraints.py``, ``_step_validation.py`` and ``_pipeline_ir.py`` are the specific evidence inspected here.

The implemented conditional IgG1 population instead conditions native draws on exact full coding RNA/protein compatibility: reject an incompatible whole draw and retry from the same parent with the next predeclared seed, at most 100 attempts. All attempts and rejection reasons are retained; no nucleotide repair, whole-J exclusion or library patch is used. This changes the unconditional native mutation distribution and is explicitly labelled conditioning, not unconstrained S5F. See ``methods guide`` and ``scripts/generate_population.py``. The generator selects validated nondefault constants explicitly through its heavy-constant and additional-reference parameters.

The tracked ``scripts/validate_constant_candidates.py`` entry point replaces the ignored exploratory helpers. It verifies every pinned raw-file SHA256 before parsing, reconstructs the IGHM genomic CDS against the transcript, and checks all four candidates’ CDS frame, annotated translation, physical protein length/hash and UniProt comparison. A fixed single coding codon tests native join arithmetic without consuming a generated variable fixture. The focused test checks all four references and rejection of tampered raw bytes. Existing raw downloads can be reused; no calculated intermediary is read, no manifest is rewritten, and no repertoire or network request is generated. Missing raw files can be restored from the manifest URLs to the listed filenames; changed upstream bytes fail instead of being silently repinned.

Reference integrity does not by itself establish population realism or assay performance.

Reusable validated records
--------------------------

``load_validated_constants(reference_dir)`` in ``scripts/validate_constant_candidates.py`` reuses the audit extraction and returns gene-keyed dictionaries. Each includes ``aa`` (the full source-specific physical constant protein), ``raw_cds_nt`` (unaltered extracted CDS), ``nt`` (the suffix to append to a native variable ending in G), ``stop_codon`` (the authentic source terminal codon), nucleotide/protein source IDs (``cds_accession``, ``protein_accession``), feature location, CDS and physical-protein SHA256 values, and the declared join convention. For IGHM/IGHG3, ``nt`` equals ``raw_cds_nt[1:]``; for IGHA1/IGHE it equals the entire extracted CDS. The overlap retains the identical native variable G rather than duplicating or replacing it.

``join_constant(variable_nt, variable_aa, constant)`` rejects incompatible terminal bases, frames and full translations; it returns exactly ``variable_nt + constant['nt']``. The audit checks this boundary with fixed ``ATGG``/``M`` inputs for all four genes. The focused test verifies exact source/hash identity, translation including the authentic stop, rejection of incompatible tails, and unchanged raw input bytes. ``validate()`` and its compact CLI output are unchanged. The existing IgG1 loader is unchanged. The population generator now selects these records explicitly through ``--heavy-constant-gene`` plus ``--additional-reference-dir``, validates before writes or native draws, and stages the selected public constant proteins for both search databases. See ``methods guide``. Fixed-input checks cover all four choices, exact complete translation and tampered-source rejection; no new native cohort, spectra or mixed-isotype frequency model was run for this integration.

Bounded native repertoire calibration
=====================================

``scripts/calibrate_repertoire.py`` executes the separately approved prospective reconstruction specification. It requires the full aggregate target TSV, explicit condition/background selection, versioned initial native model and pinned public constant-reference directories. Defaults do not launch it. Run only after source review/freeze, under the existing whole-project process/resource monitor.

The driver runs at most eight fresh 1,200-paired native generations sequentially. Related construction records are fixed at 600, with fixed seed 20261111, mutation rate 0.01, 100 prescribed read pairs per chain and equal prescribed amounts. No read files are generated. No reads, MS, Pgen, scoring or recovery selection occurs. One shared exponential NP tilt is fitted; the default fitting starts are 0, 0.10 and 0.40 inverse nt, with separately versioned starts supplied by ``--lambda-starts LOW MIDDLE HIGH``, followed by at most five widest-adjacent-interval bisections around the current best distance. Bracket widths and midpoints use exact decimal values from the stored lambda strings; equal widths choose the lower-lambda interval. Midpoints are rounded to two significant digits using decimal ROUND_HALF_UP before execution (for example, 0.175→0.18 and 0.00325→0.0033). Starts must be finite, nonnegative and strictly increasing. Adjacent intervals whose two-significant-digit midpoint repeats a tested value or reaches an endpoint are ineligible; when neither neighboring interval remains eligible, the run completes normally with ``stop_reason=parameter_precision_limit``. Otherwise it completes at eight trials with ``stop_reason=max_trials``. The receipt records all three start values, bounds and whether the selected point lies at a boundary. Every old length retains positive support. The sole objective is exact integer-histogram Wasserstein-1 distance in anchored-junction AA, retaining the full empirical target and tails. V/J discrepancies never alter the selection rule.

Target denominators must all be 33,040 genuine background records. Largest-remainder constant observation allocation preserves zero/rare bins in reporting; missing observation labels remain missing, with physical constants allocated conditionally on known calls. The allocation is deterministically shuffled independently of construction components. All full native calls remain unmasked. Full and observed-supported conditional V/J TV, residual bins, C observation/physical/imputed counts, variable/V/J length histograms, locus/family/component counts, native nonnegative integer event-count sums with valid/missing/invalid denominators (never sums of mutation rates), rejection reasons and realised SHM sums accompany every trial. Ratios retain their count denominators; source/native retry/accepted seed details remain in complete final truth.

Each native call has a 1,800-second timeout, with four hours aggregate. The outer monitor owns the combined 16-CPU, 24-GB sampled summed RSS, 64-live-thread and **10-GB whole-project** guard. Native GenAIRR generation is serial in its existing batch engine; trials are sequential because later parameter choices depend on completed results. The driver requires and reads the monitor-provided positive ``POLARS_MAX_THREADS`` allocation, clamps its positive CPU request to that ceiling and inherits the existing thread limits. It does not import scoring/vdjtools dependencies into the isolated generation environment; a missing or invalid parent allocation fails before output creation. No parallel speedup is claimed. Source, model, targets, the six used raw heavy/light constant references, declared additional-constant downloads, environment recipes and all installed GenAIRR files except bytecode/cache files (including mutation models) are pinned once and checked for drift at completion. Unresolved V categories retain their mapping-status labels; target records in zero generated bins are labelled as such, without implying catalogue absence. Generation wall time is captured immediately after the native subprocess, separately from complete trial/export time. Existing outputs are refused; no resume or cache paths exist.

Retention is declared before launch in ``execution.json``: retain deterministic gzip full native truth with verified raw/compressed hashes and byte counts; every worse trial’s model, C labels, summaries, native logs and hashes; target histograms and the complete selection record. After successful hash/roundtrip checks, remove only that trial’s newly generated template/protein FASTAs and copied reference directory. Original inputs and results are never removed. No subsequent trial consumes a previous trial’s calculated features. Failed trials remain failed; there is no replacement ninth attempt.

Prepared invocation from the analysis checkout (paths supplied by the author-owned private run specification):

.. code:: sh

   .venv/bin/python scripts/profile_build.py --output "$OUTPUT-resources.json" \
     --project-root "$PWD" --workload-root "$OUTPUT" --timeout 14400 \
     --rss-limit 24000000000 --disk-limit 10000000000 --cpus 16 -- \
     .tool-envs/generation/bin/python scripts/calibrate_repertoire.py \
     --targets "$TARGETS" --condition gel_mixture --population background \
     --initial-model "$INITIAL_MODEL" --reference-dir "$REFERENCES" \
     --additional-reference-dir "$ADDITIONAL_REFERENCES" --output "$OUTPUT" --cpus 16

.. code:: sh

   python3 tests/test_repertoire_calibration.py

The prospective v2 bracket uses ``--lambda-starts 0.4 0.8 1.6``, retaining 0.4 as a freshly rebuilt matched anchor; all trials rebuild the full native population. The earlier v1 producer and its seven valid outputs remain immutable evidence of its precision-boundary failure. This production fix does not relabel that earlier run as complete, resume it or reuse its generated features. No v2 calibration was run during implementation.

The prospective v3 numerical policy retains the same normalized NP ``(length, weight)`` pairs but supplies them in ascending ``(weight, length)`` order. This avoids small positive masses being lost when a native cumulative distribution adds them after large masses. There is no weight floor, additional rescaling, changed support or engine patch. Ordering changes inverse-CDF sampling for a fixed random seed, so v3 is a fresh realization and does not claim bitwise identity to the earlier v1/v2 anchor. One shared tilt helper owns both preflight and trial parameters. All three initial public native heavy plans must compile successfully before any output directory or native cohort is created. The receipt records compile success, numeric parameter hashes, pass structure and the ordering policy. Earlier failed roots remain untouched. V3 uses starts 0.4, 0.8 and 1.6 and a fresh ``results/fresh-matched-calibration`` output; no v3 fitting was launched during implementation.

``constant_inputs(target, seed, layout_ids=None)`` also accepts an explicit nonempty list of unique antibody IDs for a separately prepared full model. The default 1,200-record allocation and label TSV remain byte-identical to v3 inputs. Explicit 33,040 IDs use the same largest-remainder and deterministic-shuffle policy, preserving all 70 missing observation labels separately from conditional physical imputation. Construct IDs through the generator-owned ``mixture_layout``; the helper does not infer family geometry or change native annotations. Full native generation remains a separate resource-gated execution.

Existing native APIs and final outputs
--------------------------------------

- Exactly one ``annotate_truth(all_observed_chains, threads=cpus)`` call produces evaluation-only ``truth_metadata.parquet`` and ``truth_regions.parquet``, with original prescribed RNA allocations and native gene/canonical-mapping status. Native variables, constants, parent structures and truth coordinates never become observation inputs.
- Canonical observed search-entry and heavy-candidate AA sequences of at least six residues are deduplicated once. ``queries``, ``observed_aliases`` and ``excluded_queries`` retain the unique query set, every original alias/role and all rejected short/noncanonical queries. One ``exact_peptide_placements(all_truth_proteins, all_unique_queries)`` batch returns all exact I/L-compatible positions. No per-query search, N×M Python loop, worker pool or alternative matcher is added.
- One coverage level, ``arda_observed_exact_compatibility``, uses search entries only. ``observed_exact_placements`` retains the union search batch, including candidate queries explicitly decoded by aliases. ``chain_region_coverage``, ``gene_coverage``, ``locus_coverage`` reuse ``region_coverage`` and existing aggregate helpers. FR/CDR, anchored JUNCTION, separate generator CDR3, V/D/J and physical C views preserve missing/not-applicable intervals and their own denominators. Overlapping views are not summed.
- ``gene_associated_capture`` separately counts truth records and full-variable I/L classes with ANY supported protein-fragment placement, grouped by locus/segment/gene. This is not direct segment coverage or unique identity. Heavy V uses exact-mapped canonical labels only where mapping status is unique; other V and J labels remain native allele namespaces, including light alleles. C is explicitly the prescribed physical constant gene. Missing bins remain explicit. Direct templated V/J/C coverage remains in ``gene_coverage`` under its original native labels.
- ``heavy_identity_population`` and ``candidate_truth_associations`` reuse the existing full-variable compatibility criterion: a candidate placement starts at protein position zero and reaches the full variable end. Candidate/class pairs are deduplicated so physical aliases do not multiply class counts. Full-variable I/L class IDs reuse the existing ``identity_metrics`` namespace exactly: ``mass_`` plus SHA256 of the ``ms_collapse`` variable sequence. Native record IDs remain intact. Partial-variable observability requires actual overlap with variable residues: constant-only agreement is retained as an association with zero variable overlap, without being labelled partial variable capture. No dummy ranks/confirmed scores or identification recovery columns are produced. Multiple compatible classes do not establish unique source identity.
- ``catalog_status_counts`` retains all raw observed status/locus categories. Final provenance records actual candidate/search/catalog counts, query aliases, native timings, placement count, estimated residue expansions, source hashes and output hashes. It describes observability only, not MS recovery, accepted identification, experimental read weighting or a matched 1,200-record non-regression result.

Full-scale memory gate
----------------------

The reused native matcher materializes its result: the **5,000,000 placement** default guard (``--max-placements``) checks the returned count and cannot cap the engine’s internal peak before return. A monitored full-scale qualification remains necessary; there is no claimed RSS bound from toy tests or silent per-query fallback. Native annotation also materializes one whole-chain batch. The gzip raw hash is streamed, and native JSON is loaded once; raw JSON buffers are not duplicated. Native truth is released after annotation.

Coverage avoids a single whole-cohort residue explosion. It partitions deterministic whole-chain blocks (1,024 chains), computes exact placement-length plus region-length expansion estimates, and rejects any block over **2,000,000 slots** (``--max-expanded-slots``) before invoking residue expansion. Existing native Polars union/coverage runs sequentially on those coarse blocks. No chain is split, so coverage equals the unchunked algorithm; only final coverage rows are concatenated. Candidate alias-association upper bounds are checked before alias expansion as well. Limits fail explicitly rather than degrading to a hidden fallback. Profile outputs retain each block’s expansion estimate; actual wall/RSS/threads must come from the root-owned process-tree monitor. An exceeded guard needs measured review, not unrecorded expansion of the resource envelope.

Software verification and scope
-------------------------------

Constructed tests mock the annotation and placement APIs and assert a single call to each, unchanged observed-only input roles, duplicate query alias preservation, short-query exclusions, full/partial/constant-only identity semantics, empty-query behavior, coverage chunk invariance and pre-expansion rejection, stream gzip identity, amplicon and RNA canonical ``arda_fastq_native`` receipt cases, input tamper rejection, fresh-output refusal and numerically lazy module load. They do not invoke a scientific cohort or claim measured coverage/recovery. Full execution remains separately authorised and monitored; no MS or scoring stage is part of this CLI.

Simulator validation and scientific readout specification
=========================================================

Truth levels and populations
----------------------------

Maintain separate keys for family, paired antibody, chain sequence, repertoire candidate, RNA molecule/read pair, physical protein/peptide species, MS2 event, coisolated contributor and identified peptide. Heavy-chain identity is the current inference objective; light/constant molecules remain physical contributors and competitors. A paired-antibody recovery claim needs both measured pairing and an appropriate paired inference model, neither implied by heavy recovery.

Three reference populations must not be merged: (1) generated mature antibody/chain population, (2) observed sequencing-derived candidate population, and (3) actual proteomic contributors. A true protein absent from observed repertoire is not an inference false negative under a denominator restricted to observable candidates; report its absence separately and also report whole-population recovery. A true clone with no informative peptide remains in the whole-population denominator, while conditional performance is labelled explicitly.

MS2 truth is a contributor **set**, including coisolated peptides/proteins and native fractions where supplied. Do not collapse it to an arbitrary generating clone. The native identification record must link to a spectrum only after the documented count/order/RT/precursor checks pass; retain actual spectrum IDs. A correct peptide PSM is not necessarily a uniquely correct antibody assignment.

Metric dictionary
-----------------

Persist metrics as a flat long table: ``scenario_id``, ``sample_id``, ``locus``, ``region``, ``quantity``, ``value``, ``unit``, ``numerator``, ``denominator``, ``n_observations``, ``truth_level``, ``source_artifact``, ``fit_status``. Use one value per cell; unavailable numerator/denominator fields are null rather than invented zeroes. Values in count fields are exact integers. Separate chain/locus/region rows rather than packing vectors into cells. Every percentage also retains its exact numerator and denominator.

+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Quantity                              | Unit              | Numerator or definition                                              | Denominator/population                                       | Truth level                           |
+=======================================+===================+======================================================================+==============================================================+=======================================+
| Mature family size                    | descendants       | Mature descendants in one family                                     | One family                                                   | Genotype                              |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Family abundance share                | fraction          | Serum or RNA family weight, separately                               | Total weight in the same pool                                | Abundance                             |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Effective abundance diversity         | clones            | ``exp(-sum p_i log p_i)`` for normalized pool weights                | Specified RNA or serum clone pool                            | Abundance                             |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Junction length                       | amino acids       | Anchor-inclusive junction length                                     | One sequence                                                 | Genotype                              |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| IMGT CDR3 length                      | amino acids       | Anchor-excluded loop length                                          | One sequence                                                 | Genotype                              |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| SHM event burden                      | events            | Recorded mutation events along ancestral path                        | One chain/branch, specified                                  | Generative history                    |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Observable SHM burden                 | fraction          | Templated positions differing from the chosen germline               | Observed templated positions in the same chain               | Sequence evidence                     |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| V/J/C usage                           | fraction          | Chain count or pool weight for one gene                              | Same locus/pool with usable call                             | Genotype or observed candidate        |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| RNA read yield                        | read pairs        | Generated paired reads                                               | Named sample and specified input RNA count/weight            | Sequencing                            |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| ARDA retention                        | fraction          | Corrected usable candidate observations                              | Declared input/read/clone population, named explicitly       | Observed repertoire                   |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Variable-region sequence availability | fraction          | Actually observed variable amino-acid positions                      | True variable positions, for synthetic truth evaluation only | Sequencing reconstruction             |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Region proteomic coverage             | fraction          | Union of covered amino-acid positions in one region                  | Region length for that chain                                 | MS2 truth or identified PSM, separate |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Panel pooled region coverage          | fraction          | Sum of per-chain covered region positions                            | Sum of region lengths across the declared panel              | Protein evidence                      |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Peptide yield                         | peptide sequences | Distinct exported called sequences after stated normalization/filter | Named run; input weight also recorded if known               | Identified readout                    |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| MS2 yield                             | spectra           | Native MS2 event count                                               | Named run/acquisition setting                                | Instrument simulation                 |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| PSM yield                             | assignments       | Accepted PSM rows after explicit policy                              | Searched MS2 event count also recorded                       | Identification                        |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Correct peptide PSM rate              | fraction          | PSM whose normalized peptide matches a native contributor            | Evaluated accepted PSMs with truth                           | MS2 contributor set                   |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Ambiguous peptide fraction            | fraction          | Distinct observed peptides compatible with more than one candidate   | Distinct observed peptides in same candidate universe        | Sequence compatibility                |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| True heavy recovery                   | fraction          | True heavy identities meeting a predeclared accepted criterion       | All true heavy identities in the specified population        | Clone inference                       |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Conditional heavy recovery            | fraction          | Recovered true heavy identities in a named observable stratum        | True heavy identities in that stratum                        | Clone inference                       |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Identity shortlist precision          | fraction          | True heavy identities in the named shortlist                         | All candidate identities in that shortlist                   | Clone inference                       |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Target rank                           | rank              | Deterministic rank with ties recorded                                | All candidates ranked in that sample                         | Clone inference                       |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+
| Target reconstruction residual        | target units      | Synthetic summary minus recorded real summary                        | One predeclared target at matched preprocessing              | Reconstruction                        |
+---------------------------------------+-------------------+----------------------------------------------------------------------+--------------------------------------------------------------+---------------------------------------+

Define ``p_i`` as a normalized nonnegative abundance weight; use ``0 log 0 = 0``. Shannon effective diversity is descriptive, not a test. Mutation events can exceed observed substitutions because repeated changes/reversions occur; never use these burdens interchangeably. CDR3/NDN overlaps are not disjoint totals: a partition must be explicitly defined before summing region metrics. D is templated. Protein coverage is residue union, never alignment count. Report macro per-chain and pooled residue summaries separately; one does not substitute for the other.

PSM ambiguity across coisolated contributors, multiple peptide candidates and sequence-equivalent proteins needs separate labels. Do not report a native contributor intensity fraction as a posterior assignment probability. Database-search hyperscore, PSM PEP and residue ALC have different meanings; retain their native units and missingness. No new score cutoffs or multiplicity corrections are introduced here.

Realism layers and proposed parameter mechanisms
------------------------------------------------

+----------------------------------------+---------------------------------------------+-----------------------------------------------------------+------------------------------------------------+------------------------------------------------------------------+
| Parameter axis                         | Physical or biological owner                | Mechanistic intermediate                                  | Expected affected metric                       | Important distinction                                            |
+========================================+=============================================+===========================================================+================================================+==================================================================+
| Family/clone abundance concentration   | Declared repertoire scenario                | Read and protein abundance distribution                   | Effective diversity, rare-clone yield/recovery | RNA and serum pools differ                                       |
+----------------------------------------+---------------------------------------------+-----------------------------------------------------------+------------------------------------------------+------------------------------------------------------------------+
| Shared-family mutation burden/topology | Published mutation engine and declared tree | Within-family peptide sharing and discriminating residues | Ambiguity, SHM coverage, clone/family recovery | More mutations do not guarantee more usable peptides             |
+----------------------------------------+---------------------------------------------+-----------------------------------------------------------+------------------------------------------------+------------------------------------------------------------------+
| Kappa/lambda composition               | Declared paired root population             | Light-peptide species and abundance                       | Locus contribution and acquisition competition | Current tiny fixture tests kappa only                            |
+----------------------------------------+---------------------------------------------+-----------------------------------------------------------+------------------------------------------------+------------------------------------------------------------------+
| V/J/C composition                      | Reference/model and population choice       | Region sequence and cleavage opportunities                | Gene/region coverage and identity evidence     | Gene usage changes are not independent of peptide chemistry      |
+----------------------------------------+---------------------------------------------+-----------------------------------------------------------+------------------------------------------------+------------------------------------------------------------------+
| Secretion/enrichment factors           | Declared serum transfer/preparation         | Protein molar weights relative to RNA                     | Protein yield and recoverability               | Not affinity inferred from sequence                              |
+----------------------------------------+---------------------------------------------+-----------------------------------------------------------+------------------------------------------------+------------------------------------------------------------------+
| Protease identity/order                | Documented experimental digestion           | Peptide positions, lengths and cleavage products          | Junction/FR/constant coverage and identity     | Sequential digestion differs from pooling                        |
+----------------------------------------+---------------------------------------------+-----------------------------------------------------------+------------------------------------------------+------------------------------------------------------------------+
| Supported missed-cleavage behavior     | Published tool model                        | Peptide-length/species distribution                       | Yield, junction-spanning coverage              | Trypsin-trained probabilities are not calibrated for all enzymes |
+----------------------------------------+---------------------------------------------+-----------------------------------------------------------+------------------------------------------------+------------------------------------------------------------------+
| Fraction/chain recovery                | Physical preparation                        | Heavy/light/constant loading                              | Relative region yield and competition          | Loss alone cannot prove improved heavy detection                 |
+----------------------------------------+---------------------------------------------+-----------------------------------------------------------+------------------------------------------------+------------------------------------------------------------------+
| MS2 precursor capacity                 | MSSimulator native DDA selection            | Intensity-ranked competing precursors per scan/bin        | Native spectra, rare identity coverage         | Feature-only mode does not model this competition                |
+----------------------------------------+---------------------------------------------+-----------------------------------------------------------+------------------------------------------------+------------------------------------------------------------------+
| RT sampling interval                   | MSSimulator scan schedule                   | Number/timing of selection opportunities                  | Acquisition yield and repeated selection       | Seconds between scans, not Hz                                    |
+----------------------------------------+---------------------------------------------+-----------------------------------------------------------+------------------------------------------------+------------------------------------------------------------------+
| Dynamic exclusion                      | MSSimulator native selection                | Repeated versus diverse precursor selection               | Distinct peptides and clone coverage           | Parameter units and interaction must be stated                   |
+----------------------------------------+---------------------------------------------+-----------------------------------------------------------+------------------------------------------------+------------------------------------------------------------------+
| Sequencing depth/library geometry      | InSilicoSeq and actual protocol             | Observed read positions and errors                        | ARDA retention and candidate coverage          | Amplicon and fragmented RNA-seq are separate scenarios           |
+----------------------------------------+---------------------------------------------+-----------------------------------------------------------+------------------------------------------------+------------------------------------------------------------------+

These are hypotheses to inspect with the author-requested parameter variation, not guaranteed directions or measured results. Select a compact finite scenario set and parameter values before execution. Keep a declared synthetic genotype panel fixed for an assay comparison when the question is assay mechanism; regenerate the entire approved run from its pinned recipe rather than reuse computed intermediates. Repertoire-structure questions may intentionally change that panel. Do not add unrelated nulls, permutations, controls or a blind threshold/parameter scan.

Matched synthetic reconstruction contract
-----------------------------------------

1. Freeze a private target sheet from accepted retained real outputs. Record source file/section, value, units, exact counts/denominator, processing/filtering level and known missingness. Do not rerun completed analyses just to populate it.
2. Describe the real experimental design: one shared candidate repertoire versus independent subjects, reference-panel size, individual-run inventory, mixture/fraction labels and actual protocol. Preserve actual identifiers and mappings privately.
3. Specify comparable synthetic mock clonotypes by population constraints: chain lengths/regions, mutation burden, abundance/family structure and reference-panel architecture. Do not copy identifying original sequence/label artifacts into public synthetic assets by default. Reference-like does not mean sequence-identical.
4. Define which targets select parameters and which readouts remain unused during selection. Targets may include comparable filtered peptide yield, region-union coverage, protease contribution, candidate retention and accepted recovery/rank summaries where meanings match. Missing real light-chain targets remain missing.
5. Choose scientifically interpretable parameter values from actual protocol, pinned tool defaults or declared assumptions. Predeclare any optimization objective, target scaling/weights and tolerances with the author. No arbitrary numerical acceptance threshold is supplied by this document.
6. Fit only the approved targets/scenarios. Report every target residual and all attempted approved settings in the complete private record. An unused metric becomes validation only if it did not influence parameter choice; otherwise mark it fitted/post-selection.
7. Interpret reproduction as a fitted mechanistic reconstruction. Distinguish independent evidence from internal consistency; multiple settings may reproduce the same yield. Do not claim a uniquely identified protocol mechanism without evidence.

Accepted real findings and baselines remain authoritative. Definition conflicts (e.g. junction-core coverage or SHM sign) are reconciled from producing code/retained records before target use. If new work conflicts with a standing finding, tell the author before changing its presentation. No new unapproved external cohort is fetched to manufacture independence.

Figure/table specification and manual review
--------------------------------------------

Use one evidence-linked source table per metric family, with deterministic scenario and observation keys. Proposed panels: (a) locus-specific repertoire length/SHM/gene and abundance distributions, (b) observed region coverage and peptide-source composition, (c) assay parameter→intermediate→recovery response, (d) repertoire structure→ambiguity→accuracy response, and (e) real versus reconstructed target summaries plus exact settings/residual table. Panel letters here are planning placeholders, not final numbering.

A repertoire-rank plot can demonstrate abundance shape but cannot prove sequence realism; a gene-usage heatmap can show reference composition but not measurement performance. A fit-to-target plot demonstrates reproduction of those targets, not prediction. Display heavy/kappa/lambda separately, with explicit untested/missing states; do not paint absent lambda outputs as measured zero. Use native MS2 truth for identification accuracy, observed sequence masks for sequencing coverage, and known full protein only for explicitly truth-level denominators.

Acceptance and proofreading checklist
-------------------------------------

- ☐ Every requested aim has a hypothesis, measurable readout, authoritative source and prospective figure/table role.
- ☐ All biological/protocol knobs map to an actual documented tool control or explicit input composition; no homegrown surrogate model silently replaces them.
- ☐ Human heavy/kappa/lambda and secreted isotype constants are represented and separately checked; currently kappa-only smoke evidence stays labelled.
- ☐ Full RNA CDS, processed secreted protein and variable fragments have correct distinct boundaries/frame/translation; membrane tails and signal peptides cannot leak into a serum-preparation claim.
- ☐ FASTQ→ARDA→candidate coverage uses only observed sequence; missing framework is never filled from planted truth.
- ☐ Native MS2 contributor mapping passes count, order, RT and precursor checks; all coisolated contributors/fractions are retained.
- ☐ Each reported count, ratio, rank and coverage names its population/units/denominator; no PSM, molecule, peptide and clone counts are interchanged.
- ☐ Coverage is union-based and region overlaps are explicit; incomplete sequences do not shrink denominators silently.
- ☐ Truth labels/allele artifacts/source filenames are absent from scorer-visible features; candidate/readout database selection is documented.
- ☐ Accuracy separates peptide correctness, antibody ambiguity and heavy-clone recovery; tie handling and accepted criteria are preserved.
- ☐ The approved finite parameter set and target sheet precede scientific runs; no invented grid, significance scan or model selection history is omitted.
- ☐ Reconstruction marks fitted versus unused targets honestly and reports all approved target residuals without changing old real-data results.
- ☐ Fixed versions, models, seeds, resource configuration, tool provenance and source checksums identify each run. InSilicoSeq replay includes CPU count; cross-CPU FASTQ equality is not assumed.
- ☐ Every final figure/panel and rendered page is directly inspected for scientific meaning, scales, grouping, denominator and caption; tables are checked for all rows/columns.
- ☐ Final prose undergoes separate checks for scientific interpretation, numbers/provenance and readability/cross-references. No unmeasured optimality or broad first-in-field claim is introduced.

Synthetic attribution comparison
================================

Inspected inputs
----------------

All twelve original bundles have complete observed candidate columns embedded in ``{geometry}_candidate_identity.parquet``. The six leader-aware bundles also have pure candidate tables. Project the exact observed columns (sequence/NT, anchored junction coordinates, supplied genes/regions and contig/support provenance), dropping every old score, rank, mass-class/shortlist/identity label. Every original bundle retains per-enzyme final PSM tables. Their original fields correspond to ``readout._SCHEMA``, plus enzyme. Preserve original peptide/spectrum IDs, confidences and observation metadata; drop ``native_contributor_count``, ``matching_protein_count``, ``correct``, ``status`` and ``protein_ambiguous`` from scoring inputs. Never reconstruct peptide observations from only unique sequence strings.

Important correction to the preliminary specification: ``*_psm_correctness.parquet`` preserves correctness/count summaries but **does not preserve native contributor IDs or class sets**. Native truth cannot be reconstructed from those counts. Observation-level attribution evaluation needs both proteases’ retained ``ms2_truth.json`` event records (or an equivalent authoritative native contributor-ID export) for every build. These were absent from the compact local bundles and have now been selectively retrieved and verified as recorded below. Producer/native artifact receipts can locate the existing final exports for selective retrieval; no native acquisition or identification rerun is necessary or authorised.

The native truth schema used by ``simulation_metrics.psm_correctness`` is ``native_spectrum_reference`` plus ``contributors``, each containing ``peptide_aa`` and ``protein_accessions``. Truth links must be scoped by original acquisition/enzyme because scan numbers repeat between acquisitions. The original PSM ``spectrum_id`` carries that acquisition identity. Retain all coisolated contributors and all accession providers of a merged peptide. Native IDs, sequences, families and correctness labels are evaluation-only.

Proposed primary evaluation unit
--------------------------------

Use observation × full anchored heavy-junction I/L-equivalence class as the primary N-attribution unit. Physical ``generator_junction_aa`` supplies evaluation class labels; observed candidate ``cdr3_seq`` supplies candidate-class lookup without revealing that label to scoring. Distinct physical variable sequences sharing a junction are one class only in this junction-level readout, with underlying multiplicity reported. Keep existing full-variable physical/observable/recovered identities as a separate unchanged evaluation unit.

In the retained broader protein-provider view, a positive pair means the class belongs to at least one native coisolated heavy contributor at that observation. The final primary strict-N and secondary junction-interior labels are specified below. A candidate lacking a complete, concordant observed junction or compatible with unresolved class alternatives remains unresolved; partial support is not unique full-variable identity. Class scores can sum the existing spectrum responsibilities across candidate/offset explanations mapping to the same unambiguous observed junction class. Their interpretation remains conditional explanation mass, not a calibrated probability that every physical member is present. The null component and unresolved candidate mass must remain reported.

The complete physical class universe is retained per observation. Unrepresented positive classes receive zero support rather than disappearing from recall. Observations with missing native truth are unscorable, not false. Observations whose true contributors are exclusively light-chain proteins remain in false-heavy-attribution accounting; per-observation AP is undefined where there are no positive heavy classes. Pooled spectrum×junction-class AP is meaningful only if both truth labels and this class score are explicitly defined. All synthetic physical antibodies being present makes global presence-based clonotype precision an invalid discrimination benchmark.

Configurations and N eligibility
--------------------------------

Two aggregation modes (``alignment_union``, ``exclusive_spectrum``) require thirty-six fresh Evidence constructions across eighteen cases. Within each fresh build, report all contextual N minima from two through the maximum observed run plus one, retaining initial values two and three. Every configuration recreates background/search/alignment evidence and uses the original fixed synthetic ErrorModel exactly; there is no cross-configuration fitted-feature reuse. N minima are reporting/rescue eligibility parameters separate from the fixed ≥20-confirmed-nat main criterion. Peptide-identification yield and native PSM correctness are unchanged upstream observations, not claimed gains of this downstream comparison.

The existing pure ``junction_support`` module provides whole-peptide contextual run eligibility from explicit inferred nucleotide geometry. It does not call a native annotator or invent NP coordinates. Actual batched annotation must use ``annotate_junctions`` once across the complete observed-candidate junction cohort per fresh configuration, and validate returned inferred NT, native slice conventions and observed-AA correspondence. Missing/ambiguous geometry remains labelled. Whole codons must lie in one inferred N tract; mixed codons and templated D are excluded. I/L-aware agreement is required, and mismatches break a run. Existing credited observations supporting a ≥20-nat candidate are masked only from the rescue view. Germline explanations remain competitors/labels, with no empirical-null calculation.

No new birthday probability is supplied by an eligible run. The selected whole comparison family, mapped-explanation exclusions and future collision model remain distinct from these parameter readouts.

Input validation
----------------

Fixed unexplained-observation evaluation universe
-------------------------------------------------

All rescue views now include a no-N baseline (``minimum_n_run=0``, ``unexplained_rescue_no_n_filter``) and every feasible positive N minimum on the same native-evaluable observation universe minus the fresh common alignment-union explained-observation mask. Native positives are restricted to that universe before computing every rescue metric. Already resolved sources therefore cannot become rescue false negatives merely because they were deliberately masked. The separate ``all_observations_all_explanations`` gamma baseline retains the entire native-evaluable universe. Neither the complete input comparison-pair count nor the future global collision family is reduced by this rescue mask.


Explicit case selection
-----------------------

``--case-indices N...`` selects a predeclared fresh cohort and is mutually exclusive with ``--case-index``. Indices must be unique integers from 0 through 17; omission selects all eighteen. Both run provenance and execution receipts record the exact selected list. Selection never examines previous outputs and is not a resume/skip-if-exists mechanism. Every selected case still makes two fresh Evidence configurations. Original M0 qualification stays selected by unchanged ``--case-index 0``.


Fresh strict-N M0 qualification command (the output/report paths were absent at source freeze):


The raw candidate projection conditionally preserves a supplied Boolean ``background_member`` input column. It records source background-population provenance, independently of native physical truth or reference identity; it is never inferred from evaluation labels. Old scores/ranks/identity/correctness columns remain stripped. Current eighteen genuine synthetic candidate inputs do not supply this optional field, so their default scientific path is unchanged.

Fair class-score operating points, learned only on M0
-----------------------------------------------------

For a fair binary comparison, ``attribution_metrics`` now retains every complete score-tie PR point and the existing γ>0 metrics/AP, then identifies the maximum-F1 **positive class-score cutoff** on training M0 only. Equal F1 chooses the higher cutoff and thus fewer calls; ties are never broken within equal-score groups. The implicit zero-score AP endpoint is exported but is ineligible as an operating cutoff. Both the no-N baseline and every N minimum receive their own optimized score cutoff. N selection uses optimized F1 first, AP second and smaller k only to resolve exact remaining ties. Compare best N against the **optimized no-N operating point**; ranking AP and the weaker γ>0 readout remain separately visible. Main confirmed-evidence ≥20-nat calls and all producing γ scores are unchanged.

Each fresh configuration writes ``source_class_precision_recall.tsv`` with threshold, TP, FP, FN, precision, recall, F1, call count and the positive-threshold eligibility flag for every truth/view/k combination. M0 training rows retain ``max_f1_*`` fields; every case reports the applicable fixed ``operating_*`` statistics plus ``operating_threshold``, ``operating_point_available`` and ``selected_training_N``. Held-out rows do not export or promote their own maximum-F1 cutoff. All k alternatives remain retained. If a held-out case has longer feasible runs than the training grid, those extra k retain γ>0/AP/PR results but have no invented operating cutoff; every trained k is evaluated even when it produces zero calls.

The requested model parameter file is ``source_class_operating_points.json``, produced by the fresh M0 training build. Its flat metadata includes schema version, kind, exact M0 training case, strict-generative-N truth definition, fixed unexplained scope, cutoff/tie-selection rules and original source/raw-input SHA256 pins. ``operating_points`` contains one row per k, including k=0 for no-N, with the learned positive ``threshold`` (or null for no calls), training TP/FP/FN/precision/recall/F1/call count and AP. ``best_n_minimum`` names the chosen k≥2; the no-N optimized comparator is the k=0 row. No fitted Evidence/search/annotation intermediates are included. The file is an explicitly authorized learned model parameter, not a computed-feature cache.

A cohort excluding M0 must supply ``--operating-points FILE``; supplying that input while selecting M0 is rejected. The driver validates schema, training provenance, unique legal k, valid thresholds, selected-k membership and source/raw SHA pins, then includes the actual parameter-file SHA in run provenance and completion drift checks. It never fits an operating point on the other seventeen cases. Other native truth views use the same fixed training cutoffs as broader readouts. Every held-out process must consume the same M0 parameter file.


Meaningful operating-cutoff precision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Selected public operating thresholds use two significant digits in both classification and evaluation. Training considers the nearest-rounded and downward-rounded two-significant-digit boundaries derived from the original positive-score PR groups; downward rounding retains an attainable inclusion boundary when nearest rounding moves above the observed score. A single ordered sweep selects their actual maximum training F1, retaining the higher cutoff for ties. Raw full-precision PR curves and ``max_f1_*`` values remain internal training diagnostics; ``rounded_max_f1_*`` records the recomputed training counts at the deployable cutoff. The parameter JSON stores only the rounded cutoff and its corresponding counts, with ``threshold_significant_digits=2``; held-out metrics apply that exact value. Model validation rejects a supplied cutoff with excess precision. No held-out threshold optimization is introduced.

Attribution method performance comparison
=========================================

The benchmark reporter ``vdjaqp.bench.classification_metrics`` retains every candidate record and exact positive counts. Mixtures have twelve reference positives among 33,052 records; each of eight individual contexts has one held-out positive. Other records are **benchmark negatives**, not proven biological nonbinders. Supplied full-variable I/L aliases are reported separately, without silently changing the original record-label benchmark or assigning specificity to all other sequences.

Primary metric is average precision (AP): sum over distinct descending score thresholds of recall increment times precision at that threshold, with tied scores entered together. This is not trapezoidal PR area. ROC-AUC integrates the grouped ROC curve and gives tied positive/negative pairs half credit. Deterministic identifier tie-breaking is used only for reported ranks/top50/top500, not AP/AUC. The exact floating median is preserved. Precision, recall and F1 use the existing ``confirmed_nats >= 20`` call; no optimisation or threshold scan is performed. Empty shortlists have precision zero. With no benchmark negatives ROC-AUC is unavailable.

The fresh driver compares two accepted requested final baselines (original and corrected arithmetic, both alignment_union) with corrected exclusive_spectrum. Corrected alignment_union is the causal aggregation comparator; original arithmetic remains the fixed original acceptance baseline. All ten contexts and every tradeoff remain visible. It does not reuse alignments/error models between fresh variants. Mixtures use existing individual-run calibration; each individual run is split by its original source stem and excludes its own target from calibration.

Bounded native assay-parameter comparison
=========================================

Two explicit published parameters
---------------------------------

``scripts/smoke_proteomics.py`` accepts:

+---------------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------+-----------------+
| Option                          | Native MSSimulator parameter                                                                           | Domain / meaning                                                                                               | Default         |
+=================================+========================================================================================================+================================================================================================================+=================+
| ``--min-detectability``         | ``algorithm:MSSim:Detectability:min_detect``                                                           | Finite score cutoff in [0,1]; retain predicted peptide detectability **strictly greater** than this cutoff     | 0.5             |
+---------------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------+-----------------+
| ``--dynamic-exclusion-seconds`` | ``algorithm:MSSim:RawTandemSignal:Precursor:Exclusion:use_dynamic_exclusion`` and ``…:exclusion_time`` | Finite nonnegative seconds; zero sets exclusion false; a positive duration enables it and sets the native time | 0, disabled     |
+---------------------------------+--------------------------------------------------------------------------------------------------------+----------------------------------------------------------------------------------------------------------------+-----------------+

The native inactive ``exclusion_time`` remains unchanged when the CLI duration is zero; the existing native default is 100 seconds. Setting zero does not rewrite that inactive native field to zero. Missing native keys fail the existing complete-key assertion instead of silently ignoring an option. The CLI rejects invalid values before native execution or output creation; ``forward_draw`` validates them before reading inputs or creating output. The forward CLI propagates both options unchanged to all three acquisitions and records them in ``discrepancy.json``; ``joint_draw.json`` records them with all native commands and per-acquisition summaries. The default settings preserve the existing simulation; the extra provenance fields are explicit additions.

Published source pin: OpenMS 3.0.0, ``60ce8966d577b43a6419403dec171539e08fe513``. `DetectabilitySimulation.cpp:194–240 <https://github.com/OpenMS/OpenMS/blob/60ce8966d577b43a6419403dec171539e08fe513/src/openms/source/SIMULATION/DetectabilitySimulation.cpp#L194>`__ predicts from peptide sequence, applies the strict cutoff, and declares the 0.5 default. `OfflinePrecursorIonSelection.cpp:61 <https://github.com/OpenMS/OpenMS/blob/60ce8966d577b43a6419403dec171539e08fe513/src/openms/source/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.cpp#L61>`__ defines the exclusion duration in seconds. `Lines 433–440 <https://github.com/OpenMS/OpenMS/blob/60ce8966d577b43a6419403dec171539e08fe513/src/openms/source/ANALYSIS/TARGETED/OfflinePrecursorIonSelection.cpp#L433>`__ convert duration to integer scan intervals using floor(duration/RT interval); disabled exclusion only prevents repeated selection within one scan.

The 0.75 condition is a **learned-score stringency scenario**, not a measured physical sensitivity or calibrated identification probability. Thirty seconds follows one available acquisition-method description; the instrument applicable to the experiment remains unresolved, so it is a **reported-method scenario assumption**, not a verified instrument setting. No other settings from that description are silently imported. In particular, the simulator’s five selections per RT bin are not equated to the experimental instrument’s top-N/cycle time.

Finite design and selection rule
--------------------------------

+-----------------------+----------------------+--------------------+---------------------------------------------------------+
| Label                 | Detectability cutoff | Exclusion duration | Execution                                               |
+=======================+======================+====================+=========================================================+
| baseline              | 0.5                  | off                | Accepted saved equal-amount qualification; do not rerun |
+-----------------------+----------------------+--------------------+---------------------------------------------------------+
| exclusion30           | 0.5                  | 30 s               | One fresh joint draw                                    |
+-----------------------+----------------------+--------------------+---------------------------------------------------------+
| detect075             | 0.75                 | off                | One fresh joint draw                                    |
+-----------------------+----------------------+--------------------+---------------------------------------------------------+
| detect075_exclusion30 | 0.75                 | 30 s               | One fresh joint draw                                    |
+-----------------------+----------------------+--------------------+---------------------------------------------------------+

Every new joint draw uses twelve relative amounts of exactly one, summing to twelve; 10,000 native intensity units per physical record and 120,000 total; the same supplied twelve physical protein records and all 33,052 search candidates; the same 1,200-s gradient, two-second sampling, five MS2 selections/bin, published detectability model, fixed theoretical tandem mode, search settings, RNG settings, software/environment and peptide definitions. No sequence completion, substrate replacement, noise calibration, new seed, abundance refit, extra setting or control is introduced. Source/raw/executable hashes and the accepted qualification hashes must accompany the final comparison.

The existing joint API produces solution trypsin, gel trypsin and one shared chymotrypsin acquisition, regenerating every block for each new draw. This gives exactly **three new joint draws / nine fresh acquisitions**, not twelve new acquisitions. The two identically configured trypsin acquisitions may predict identical values under the fixed reproducible native RNG; they are not modeled as independent stochastic experimental realizations. Their distinct observed targets/residuals remain. Shared chymotrypsin enters the selection vector once. No preparation-specific chemistry is invented.

Select by the **existing 255-coordinate discrepancy only**: three compatible-species yields and 3×12×7 per-reference region-coverage fractions, with the existing log-yield transformation, weights and denominators (``vdjaqp/abundance_inference.py:target_vector``, ``distance``). Rank all four conditions by that scalar distance; retain every per-coordinate signed residual and every setting, including worse alternatives. For an exact distance tie, prefer the accepted baseline, then the order ``exclusion30``, ``detect075``, ``detect075_exclusion30``. This tie rule is bookkeeping, not evidence of a measured advantage. Final discrepancy is an in-sample reconstruction criterion and cannot establish absolute experimental calibration or optimality outside these four settings.

Keep recovery/identification accuracy, ≥20-confirmed-nat calls and source ambiguity **out of parameter selection**; report them afterward as evaluation readouts without changing thresholds. Report all endpoints and denominators for species, per-reference residue coverage, PSM truth linkage and candidate/class recovery. Lower discrepancy alone does not guarantee better identification or a physically validated assay model. The observed versus native search/export observation-process difference remains; no arbitrary native-intensity-to-experimental-Area threshold is added.

Concrete commands after source review
-------------------------------------

Run from the repository root. The invoking private analysis supplies absolute ``ASSAY_REPERTOIRE``, ``ASSAY_TARGET_IDS``, ``ASSAY_TARGETS``, and ``ASSAY_SOURCE_FORMAT`` from its approved input/source records; do not put stakeholder paths or identifiers in this public-ready document. Use a new ``ASSAY_OUTPUT`` subtree; existing outputs are never resumed or overwritten. The commands below run only the three specified new conditions. The accepted baseline qualification remains an input deliverable for the final comparison, never a feature cache or an execution input.

Matched-reference MS bridge
===========================

``scripts/run_matched_reference_ms.py`` provides three explicit gated stages: ``qualify-regions``, ``acquire``, and ``evaluate``. It does not generate a repertoire/reads, resume work, fit parameters, alter core scoring, or run another stage automatically. Each stage is invoked explicitly. Source/query packets remain private local inputs; no private identifying labels or fixed stakeholder paths are embedded in this public-ready CLI.

Approved preparation inventory
------------------------------

Each of the two geometries has **three unique acquisitions**, six acquisitions overall:

=================== ================= ================= ============
Block               Physical H chains Physical L chains Enzyme
=================== ================= ================= ============
solution_trypsin    12                12                Trypsin
heavy_band_trypsin  12                0                 Trypsin
shared_chymotrypsin 12                0                 Chymotrypsin
=================== ================= ================= ============

Per geometry, solution trypsin plus the one shared chymotrypsin produces the solution pool; heavy-band trypsin plus that same shared acquisition produces the heavy-band pool. Thus there are **four pooled scores**, with shared chymotrypsin entering each pool once. No acquisition crosses geometry boundaries. This mirrors the already approved heavy-band/shared-file experiment, not a new gel-benefit contrast. Heavy intensity stays 10,000 per represented chain in every preparation; unchanged absolute heavy-band recovery is an explicit assumption, not a measured recovery/enrichment factor. Light chains in heavy-band outputs are unassayed, not measured zeros.

Fixed baseline is 1,200 s, five MS2/bin, learned detectability 0.50, dynamic exclusion off, equal input amounts and conventional equal H:L represented-chain ratio in solution. The CLI exposes no sweep, abundance fitting or recovery-driven parameter choice. The qualified 24-chain packet is unchanged; heavy-only raw model views are derived from its exact H rows, without altering variable/constant sequence or constructing new partners.

Source-region stage
-------------------

.. code:: sh

   .venv/bin/python scripts/run_matched_reference_ms.py qualify-regions \
     --packet-manifest /private/qualified_packet/packet_manifest.json \
     --output /fresh/reference_regions --cpus 16

The driver verifies packet artifact/source hashes and makes **one** published ARDA AA annotation call on all 24 reported variable proteins, ``reannotate_contigs(..., seqtype='aa', map_d=False)``. It never reverse-translates or invents GenAIRR NT/genotypes. The installed ``arda.annotate.transfer`` AA branch defines region slices as ``query[qs-1:qe]``, establishing one-based inclusive AA inputs. The adapter converts those bounds to zero-based half-open intervals and checks source slices/ranges. Native best-call ambiguities remain native calls, not canonical-gene proof. Raw annotations are final retained outputs. Supplied V/J labels remain separate metadata from inferred V/J fields and never overwrite heavy query fields.

Physical full-variable and C intervals follow the unchanged reported variable/public-constant join. Heavy anchored junction comes from the preserved source query; light junction uses a uniquely placed native-inferred AA substring when available. Junction includes Cys and terminal F/W; CDR3 excludes anchors. Inferred FR3/CDR3/FR4 endpoints conflicting with a known junction become explicitly unavailable rather than repaired. D, NP1/NP2 and genealogy are unavailable for this reported AA-only source, not borrowed from generated truth. ``physical_metadata.parquet``, ``physical_regions.parquet``, raw AA annotations and a source-pinned manifest are private final evaluation/source deliverables.

Acquisition interface
---------------------

.. code:: sh

   .venv/bin/python scripts/run_matched_reference_ms.py acquire \
     --packet-manifest /private/qualified_packet/packet_manifest.json \
     --reference-regions /final/reference_regions/manifest.json \
     --observed-receipt /final/amplicon/observed/receipt.json \
     --geometry amplicon --block solution_trypsin \
     --mssimulator /verified/runtime/MSSimulator --stage-timeout 900 \
     --output /fresh/amplicon_solution_trypsin --cpus 16

The canonical native ARDA ``observed/receipt.json`` pins the three exported frames, actual counts, source/code hashes, locus counts and background roles. Direct FASTQ input uses the native ARDA route only; there is no separate mapper/contig support receipt. The observed-row provenance distinguishes uncorrected reads, native assembled reads and native-corrected junctions without replacing their native NT/AA values or changing default quality settings. The append-only API preserves all actual observed background rows and adds exactly twelve supplied queries with ``background_member=False``. It preserves supplied AA, reconstructed NT, regions and calls; the separate scoring view alone uses existing gene normalization. The actual heavy candidate denominator is N + 12, never a forced 33,052.

Search FASTA contains only support-eligible observed entries, exact supplied H query sequences, reported VL-only auxiliary entries and independently declared public constant competitors, with identity-collision rejection. The full physical proteins are NEVER a candidate database or inferred candidate-boundary source. Physical simulation and native ``smoke_proteomics.py`` use the declared per-block raw model view, explicit fixed settings and pinned MSSimulator. Auto-fetch is disabled. Native output/linkage is validated with existing ``read_idxml`` and ``psm_correctness``; no native hyperscore is converted into confidence probabilities.

After complete idXML/TextExporter/native-contributor and parsed-call validation, the driver hashes and retires only its own ``spectra.mzML`` and consumed ``truth.featureXML``. The ledger records raw SHA/bytes, validation time and logical task bytes before/after; whole-project allocated-disk measurements belong to the outer resource receipt. Producer summary files remain unchanged. Native final idXML, truth idXML, TextExporter TSV and contributor JSON are deterministically gzip-retained with raw/compressed round-trip SHA/bytes. Raw-call/correctness/accuracy Parquets are final readouts, not fitted feature caches. Failures retain diagnostics and do not trigger deletion/retry. Source FASTAs/query frames and requested native final deliverables are retained.

Pooled evaluation stage
-----------------------

.. code:: sh

   .venv/bin/python scripts/run_matched_reference_ms.py evaluate \
     --packet-manifest /private/qualified_packet/packet_manifest.json \
     --reference-regions /final/reference_regions/manifest.json \
     --observed-receipt /final/amplicon/observed/receipt.json --geometry amplicon \
     --acquisition-manifests /final/solution_trypsin/manifest.json \
       /final/heavy_band_trypsin/manifest.json /final/shared_chymotrypsin/manifest.json \
     --background-truth-gzip /private/native/truth.json.gz \
     --background-validation-receipt /private/validation/summary.json \
     --output /fresh/amplicon_pooled_evaluation --cpus 16

The three final acquisitions must have unique declared block identities, the same geometry, matching observation/packet/region pins and correct per-block physical H/L counts. Raw idXML is reparsed with original run IDs; no old alignment, calibration or scored rows are loaded. Two complete fresh ``Evidence.from_frames`` builds use the unchanged simulation ErrorModel (database error 0.10), existing aggregation/isotype prior and ≥20-confirmed-nat rule. Background truth and physical AA annotations are never arguments to inference.

A single exact evaluation batch combines all unique candidate AA, raw called-peptide AA and native contributor AA of at least six residues, preserving excluded short strings. It searches all raw generated background H proteins plus the 24 physical reference proteins. No computed observation features substitute for the raw background truth; its compressed/raw identity, complete native validation, model, observed-chain and heavy counts are checked. Candidate aliases are deduplicated before search and restored only after bounded expansion estimates. Returned native placements are evaluation outputs, not inference input.

Optional matched experimental yield/coverage report
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Add ``--experimental-targets /private/recorded/experimental_fit_targets.tsv`` to ``evaluate`` to produce the requested experimental comparison. This optional input is source-hashed and rechecked with the other stage inputs; old output directories remain immutable. It never enters candidate construction, native acquisition, the scorer, recovery calls or parameter selection.

The report reuses existing exact full-protein placements, remaps I/L-equivalent query aliases to one species, and retains only whole peptides ending within the supplied heavy variable sequence. C-crossing peptides are excluded from this comparison rather than clipped. The twelve unchanged supplied-query masks define matching FR1–FR4/CDR1–CDR3, anchored junction and variable regions; native physical-AA-inferred masks remain a separately labelled evaluation. Missing or mismatched reference keys/residue denominators fail explicitly.

For each of the three unique acquisitions and two condition pools, ``*_query_species.parquet`` retains native MS2 group counts, PSM assignment rows, noncanonical/short diagnostic flags, distinct eligible species, variable-panel compatible species and unmapped species. Noncanonical and short flags can overlap and must not be summed as sequential losses. Pooled species and coverage are unions; native group counts sum the two unique within-pool acquisitions. ``*_supplied_query_coverage*.parquet`` retains all twelve-reference/nine-region coverage denominators and pooled residue unions. These compatibility measures do not establish unique peptide source or spectrum identification accuracy. The original experimental PEAKS exports have their documented upstream Area/ppm/PTM processing; no analogous Area, ALC or hyperscore cutoff is invented for native search calls.

``experimental_target_residuals.tsv``/``.parquet`` retain per-reference and pooled signed/absolute residuals separately for solution trypsin, heavy-band trypsin, shared chymotrypsin and both condition unions. The shared chymotrypsin display aliases retain identical predictions. ``experimental_comparison.json`` uses the existing 255-coordinate target and distance definitions: three species counts plus 3 × 12 × 7 supplied-region fractions, with shared chymotrypsin counted once per geometry. It retains predicted and recorded vectors, target hash and ``recovery_used_for_selection=false``. Across both geometries the inventory remains six acquisitions and four pooled scores; no acquisition, exact search, scoring build, fit or scientific contrast is added by this report.

Distinct accuracy units
~~~~~~~~~~~~~~~~~~~~~~~

1. **Canonical query ROW retrieval** uses the existing ``classification_metrics`` directly on every candidate row, twelve original supplied query IDs as positives. Other rows are explicit benchmark negatives; this endpoint retains the original-compatible label definition and does not imply absence of an AA-equivalent physical class.
2. **Physical H full-variable I/L class recovery** reuses existing ``identity_metrics`` on the twelve physical H chains and exact evaluation associations. Partial fragments/shared ambiguous full matches do not establish unique recovered physical identity. Physical rows, distinct classes and database-observable classes retain distinct denominators.
3. **End-to-end CLASS presence** uses the union of all raw-generated H variable I/L classes and known physical H reference classes, with identical classes collapsed. Positive classes are physical references only, including any raw-background alias; other raw-background classes are absent from the declared physical assay. Full compatibility requires candidate placement start zero/end at or beyond the actual variable end. A candidate credits a class only when compatible with exactly one full class. Existing max score and max confirmed nats are aggregated independently per class. Ambiguous and no-full candidates remain unresolved, not silently assigned.

Missing/unobservable classes receive an explicit evaluation-only ranking sentinel ``min(observed finite class scores)-1``, or -1 if no class is observed, and zero confirmed nats. This places missing classes below all observed negative as well as positive finite scores; it never changes candidate inference, thresholds or native parameters. Exact class ID is the existing ``mass_`` + SHA256 of I/L-collapsed full variable AA. Existing classification metrics compute P/R/F1/AP/tie-half-credit ROC-AUC on class-ID rows. Report both the complete end-to-end universe and the separately labelled full-compatibility database-observable subset; ambiguity may make a database-observable class uncreditable. All class/population counts are retained. No metric is promoted as clone posterior calibration.

A separate **all-shortlist candidate-record precision** denominator includes every called candidate, including unresolved ones, and reports definitive positive, definitive negative and unresolved called-record counts. It does not mix candidate counts with class counts or hide unresolved calls in a class-only precision denominator. Full class-association and candidate-class-mapping tables remain final deliverables.

Coverage and memory
~~~~~~~~~~~~~~~~~~~

Per-condition physical AA-region coverage distinguishes native contributor compatibility, called-peptide compatibility and accepted-heavy-candidate exact compatibility. The latter is compatible sequence support, not L/paired identification. V/J/FR/CDR are source inferred; full C is the declared model. Heavy-band L rows have null coverage/denominators and ``not_assayed`` status; raw source intervals remain in the separate region packet. Physical-assay inventory states the 24/12/12 chain memberships explicitly. Existing ``capture_tables`` retains all/positive-LLR/accepted candidate gene/junction proxy capture with its original interval definitions, without experimental read weighting.

The default native-placement and alias-association bound is 5,000,000 rows. The existing engine materializes hits before that returned-count guard, so monitored native peak qualification remains mandatory; no hidden per-query fallback or performance claim follows from tests. Existing bounded coverage uses whole-chain blocks and a 2,000,000-expanded-slot bound before residue expansion. All native searches remain batched; no N×M Python matching or worker pool is introduced.

Missing isotype evidence
------------------------

Generic producer calls ``IGHG`` and subclass calls ``IGHG1`` through ``IGHG4``
are IgG. Missing C calls and unresolved ``IGHC`` remain ``unknown``.
When unknown calls occur, they receive zero isotype log-evidence. Known
categories receive the logarithm of target-category probability divided by
observed-category frequency, with frequencies computed over known calls only.
This assumes the known calls describe the available category mix; it does not
establish that missingness is random. A single observed category supplies no
relative category information. Fully annotated inputs retain the original
score calculation. Variant candidates retain native corrected-clone C-call
provenance separately from the observation carrying their protein sequence.

Experimental conditional source-class presence
-----------------------------------------------

``vdjaqp.aggregate.compete_identity_classes`` and ``fit_identity_mixture`` are
optional public research APIs. They do not replace ``Evidence.score`` or change
the original confirmed-evidence acceptance criterion. Their source-presence
operating point requires independent training and held-out validation; the
implementation alone establishes neither improved accuracy nor superiority to
another method.

The caller supplies an explicit mapping from every scored candidate to an
identity class. For an anchored-junction endpoint, conserved anchors and I/L
equivalence must be resolved consistently by the caller before constructing this
mapping. The APIs do not infer receptor boundaries, physical membership or truth
labels. Candidate aliases and alternative alignment offsets within one class do
not create additional class likelihood terms.

For spectrum :math:`s` and class :math:`c`, the class log weight is the maximum
across its compatible scored placements :math:`a`:

.. math::

   \ell_{sc}=\max_{a\in(s,c)}
   \{\mathrm{LLR}_a+2\log q_a-\log r_a-\log n_a\}.

Here :math:`\mathrm{LLR}_a` is the existing alignment log likelihood ratio,
including its unchanged substitution/error and mass treatment; :math:`q_a` is
peptide confidence, :math:`r_a` the alternative-call rank and :math:`n_a` its
read multiplicity. Confidence follows the existing responsibility convention:
missing/NaN confidence becomes 0.5 and confidence is clipped to [0.001, 1]; rank
and read multiplicity are floored at one. A class with no compatible placement
has likelihood zero for that spectrum. Taking a maximum is an explicit
best-placement approximation, rather than a sum over aliases or offsets.
``compete_identity_classes`` also reports equal-weight class-versus-unit-null
responsibilities and a unique-MAP diagnostic. These are not the fitted mixture
weights or calibrated probabilities of source presence.

The mixture includes a unit-likelihood null component and nonnegative weights
:math:`w_0,w_1,\ldots,w_C` summing to one. It fits observed likelihoods only:

.. math::

   G_s(w)=w_0+\sum_{c=1}^{C}w_c\exp(\ell_{sc}),
   \qquad Q(w)=\sum_s a_s\log G_s(w).

For each spectrum, the sorted set of distinct I/L-equivalent alternative peptide
sequences defines its observation species. If that same set occurs in
:math:`N_s` spectra, :math:`a_s=1/N_s`; its repeated acquisition therefore has
total unit fitting weight. Different alternative sets remain different species.
This observation-weighting convention does not estimate protein abundance.
Weights start uniformly across matched classes and the null, and are updated by
expectation maximization. The default stopping rule is maximum absolute weight
change at most :math:`10^{-8}`, with at most 1,000 iterations. The convergence
receipt retains the objective history and numerical status. Sparse native-array
reductions avoid an observations-by-classes matrix. The null weight is also
retained in log space to avoid losing its contribution after linear underflow.
Multiple classes can retain positive weights; no per-spectrum winner is required.
The model is a conditional mixture of class explanations, not an explicit
physical model of simultaneous fragment production by coisolated precursors.

For each fitted class, its removal contrast holds the other fitted weights fixed
and renormalizes them:

.. math::

   G_{s,-c}=\frac{G_s-w_c\exp(\ell_{sc})}{1-w_c},
   \qquad \Delta_c=\sum_s a_s\{\log G_s-\log G_{s,-c}\}.

All observations enter this contrast, including those without a placement in
class :math:`c`. Stable log sums avoid subtracting a dominant contribution from
a nearly equal full denominator. Other weights are **not** refitted after
removal: :math:`\Delta_c` is neither a profiled likelihood-ratio statistic nor a
calibrated posterior, and has no automatic significance interpretation. Near
boundary solutions can depend on the fixed fitting/stopping convention.
Unavailable contrasts must not become acceptance calls. Nonconvergence remains
explicit; a fixed-other contrast from such a fit does not certify presence.

Classes with exactly equal sparse likelihood profiles across all observed
spectra form an observational equivalence group. Equality uses the same spectrum
keys and finite input log weights, with missing likelihoods consistently zero;
no additional numerical tolerance is introduced. The returned
``identity_equivalence_group``, ``identity_equivalence_size`` and ``identifiable``
fields retain this ambiguity. An indistinguishable group cannot establish its
individual members merely because their fitted weights or removal contrasts are
positive. A unique profile is a necessary identifiability check, not proof that
the candidate is physically present. Group identifiers are local to one fit. Distinct individual profiles can still
have nonidentifiable joint explanations: separate V1/V2 and J1/J2 peptide
witnesses may support V1J1 + V2J2 and V1J2 + V2J1 equally. A fixed-other removal
contrast can be positive in this situation because it does not refit the other
classes. The exact-profile guard does not resolve this global mixture ambiguity;
a positive contrast cannot by itself establish an individual member of these
alternative explanations.

The fit returns matched classes. Callers must left join the complete declared
candidate-class catalog, assigning zero weight and contrast only to genuinely
unmatched classes, and retain unresolved groups separately. A class-presence
benchmark must define its shared identity unit, preserve all physical positives
including unobserved sources, and distinguish physical from observable recall.
Any contrast threshold is selected on training simulations only and frozen before
held-out evaluation, together with the software, error model, weighting and
convergence recipe. Full source populations, adverse counts, ambiguity and the
original scorer's baseline remain reportable. A small number of simulated
physical panels limits generalization; repeated spectra or classes within a
panel are not independent validation panels.

Certified conservative profiled-removal bounds
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The optional ``profile_identity_removals`` API refits the remaining mixture
weights after each requested class removal, holding the observed likelihood
matrix, background/error model and observation weights unchanged. It provides a
conservative bound rather than treating an unfinished numerical fit as an exact
profile optimum. The original scorer and fixed-other diagnostic remain intact.

Let :math:`A=\sum_s a_s`, :math:`g=\nabla Q(w)` and
:math:`U(w)=\max(0,\max_k g_k-A)` at any feasible fitted weight vector. Since
:math:`g\cdot w=A` and :math:`Q` is concave on the weight simplex,
:math:`Q(w)\leq Q^*\leq Q(w)+U(w)`. Sparse reductions compute every gradient
component, including the null, without a dense likelihood matrix. Let
:math:`\widehat Q_F` and :math:`\widehat Q_{-c}` be the full and reduced fitted
objectives and :math:`U_F,U_{-c}` their optimization-gap bounds. Then

.. math::

   B_c=\max(0,\widehat Q_F-\widehat Q_{-c}-U_{-c})
   \;\leq\; Q_F^*-Q_{-c}^* .

Neither fitted point must be an optimum for this lower bound to hold. The
certificate-based rule therefore retains convergence status, objective and
gradient gaps explicitly rather than requiring convergence as a surrogate for
validity. Nonfinite numerical certificates and failed required refits remain
unavailable. The fixed 1,000-step/weight-change fitting recipe is recorded before
training calibration; it is not tuned against held-out truth. A loose gap can
reduce sensitivity and does not establish source absence.

The exact profiled loss is bounded above by
:math:`\Delta_c+U_F`. The conservative score :math:`B_c` is additionally bounded
above by :math:`\max(0,\Delta_c)`: the reduced optimum cannot be worse than the
feasible fixed-other removal. These bounds can exclude a class from exceeding
a predeclared score threshold without an expensive reduced refit. Screened
classes retain the bound and status; unprofiled classes without a sufficient bound
remain unresolved, never silently assigned an accepted zero/error surrogate.

Acceptance in this experimental variant uses an available finite :math:`B_c`
strictly above a training-selected, frozen nonnegative threshold. Exact equal
profiles remain individually unresolved. The four-class unphased V/J example
has positive fixed-other contrasts but zero certified lower bounds because
alternative mixtures explain the observations. These are conditional model
evidence bounds, not calibrated posterior probabilities, significance tests,
protein-abundance estimates or proof of universal biological identifiability.

Optional observed-AA emission with a common class-average null
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``AlignmentScorer.score(..., include_components=True)`` and
``Evidence.from_frames(..., include_components=True)`` retain three optional
alignment components in the same scoring call. ``log_context_emission`` and
``log_context_null`` decompose the original context-masked information LLR.
``log_emission`` instead uses observed amino-acid compatibility and the unchanged
error/substitution model, independently of germline-context availability. A known
amino-acid mismatch therefore remains a mismatch even when its gene coordinate
is unavailable. Noncanonical aligned observed/source residues make this emission
unavailable (NaN), rather than supplying an artificially favorable incomplete
likelihood. The original ``llr`` values, default output schema and default
aggregate scores are unchanged.

The original information LLR uses a candidate-position-specific background.
Exponentiating it does not make comparable class likelihoods relative to one
shared null: a rarer candidate context can outweigh closer observed sequence
compatibility. The experimental ``compete_identity_emissions`` API instead starts
from the observed-AA emission components. For peptide :math:`p` and class
:math:`c`, let :math:`E_{pc}` be the largest emission among its supported
placements and let :math:`C` be the full declared class count. The common null is

.. math::

   N_p=\frac{1}{C}\sum_{c=1}^{C}E_{pc}.

Missing search-supported placements contribute zero to this average; candidate
aliases and alternative offsets do not multiply a class term. The API constructs
class ratios :math:`E_{pc}/N_p` and then applies the existing confidence/rank/read
convention before spectrum-level class competition and the optional certified
removal calculation. It returns per-peptide null-construction metadata alongside
the competition table. Peptides without any supported class are recorded as
null-only and do not create source-class evidence.

This uniform class-average emission is an explicit candidate-pool background
approximation. It is conditioned on the declared search support and is not a
calibrated distribution of incorrect database-search assignments, an observation
of protein abundance or a replacement for search-score reliability calibration.
Physical source labels do not enter its construction. Raw emissions, original
information LLRs, null metadata and all declared candidate classes remain
separately reportable. Validation must use freshly produced observations and a
training-only operating point, with the complete software-component recipe
frozen before held-out evaluation.

Optional native-root-consistent observed representatives
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``scripts/smoke_repertoire.py --representative-policy native_root`` and
``export_candidates(..., representative_policy='native_root')`` optionally prefer
an already observed complete frame whose observed junction nucleotide sequence
exactly matches its corrected native parent's junction. Length, annotation score
and the existing stable read/row/assembled tie-breakers then resolve alternatives.
The lazy grouped selection joins only native membership and corrected sequence;
protein identities, peptide observations and physical-source labels do not enter.

The default ``longest`` policy preserves the original selection and output
schema. The optional policy never replaces bases with the corrected junction,
constructs a consensus, fills missing sequence or filters parents by read count.
If no complete root-matching frame exists, the original representative selection
is retained. Corrected-junction-only fallback records remain distinct from
observed frames, and all original native variants remain in their source files.

Optional exports add ``representative_policy``,
``representative_native_root_match``, ``native_root_matching_observations`` and
``representative_selection_status`` to the catalog, candidates and search
entries. A root-matching observed frame is labelled ``observed_native_root``;
a complete observed fallback is ``observed_fallback_no_root_match``. A
``junction_only_no_complete_observation`` record has a null observed-root-match
flag rather than claiming observed sequence support. The receipt records the
policy and selection-status counts. Existing parent/background identities and
native observation counts retain their original meanings. Search-entry eligibility
retains the original parent's maximum complete-frame length, or the corrected
junction length when no complete frame exists. A preferred short observed frame
therefore remains in its original search entry without invented extension;
peptide matching cannot supply missing extent. Parent read support
is not support for every selected full-sequence variant or a protein-amount
measurement. Selection is a distinct preprocessing choice requiring fresh
matched evaluation, not a change to the original protein evidence criterion.

Optional typed all-locus candidate export
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``scripts/smoke_repertoire.py --all-loci`` and
``export_candidates(..., all_loci=True)`` additionally write
``all_locus_candidates.parquet`` with the existing typed candidate fields for
each usable native IGH, IGK and IGL parent. The heavy-only candidates, catalog
and search entries retain their original values and schemas. The optional table
uses the same observed-frame selection and native membership as those outputs;
the receipt pins its hash and row count. It does not extend observed sequence,
infer chain pairing or turn junction-only records into complete variable domains.
The exporter retains its requirement for a usable heavy candidate. Downstream
all-chain inference and per-chain calibration require an explicit matched design;
this export alone does not validate transport of heavy-chain calibration.

Optional RNA-junction reliability penalty
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``identity_rna_support(candidates, identities)`` requires the complete observed
``native_root`` candidate population and its anchored I/L-equivalent junction
class mapping. Aliases of a native parent must agree on class and native count.
Class support is the maximum qualifying root-matching observation count over
native parents, never a sum over aliases or parents. These raw/assembled
observations can overlap and are not independent molecules, whole-sequence
quality estimates or protein amounts. Missing metadata fails closed; zero-support
classes remain in the population.

``fit_rna_junction_reliability(support, genuine)`` optionally fits junction
reliability from one row per observed training class, labelled against *all*
generating RNA classes rather than the physically selected protein subset. At
each tied support, the error estimate is ``(errors + 0.5)/(classes + 1)``.
Repeating this estimate once per original class before the existing
non-increasing isotonic fit preserves class-count weights and identical fitted
values within ties. The function returns unique support knots and genuine-class
probabilities strictly between zero and one. Prediction interpolates linearly
between knots and remains constant beyond their endpoints.

``apply_rna_junction_reliability(class_scores, support, knots, q)`` takes only
frozen calibration and observed support, with no truth argument. For the existing
nonnegative certified removal bound ``B`` and interpolated reliability ``q``, it
appends ``S = max(0, B + log(q))``. Since ``0 < q <= 1``, ``0 <= S <= B``;
``q = 1`` reproduces the original bound exactly. The helper preserves original
scores, class rows, ambiguity flags and unavailable-score states. Its added
``rna_adjusted_score_available`` flag never makes an unavailable bound available.
The penalty is neither a protein-presence posterior nor an FDR estimate; its
calibration is conditional on the training RNA depth and native correction model.
An operating threshold must be calibrated on training panels and frozen before
independent validation, retaining the unpenalized criterion as a matched baseline.

Conditional repertoire sampling display
---------------------------------------

``scripts/simulate_repertoire_sampling.py`` generates fresh multinomial draws
from explicitly fitted per-library junction frequencies at observed UMI-molecule
depths. It includes an unresolved remainder category and retains shared identity
keys across libraries. Its intervals describe conditional sampling variability;
it does not infer unseen clones or simulate RNA/MS identification performance.

.. code-block:: sh

   .venv/bin/python scripts/simulate_repertoire_sampling.py \
     --abundance /private/retained/junction_abundance.parquet \
     --processing /private/retained/processing_plot.tsv \
     --output /private/new_sampling_output --replicates 200 --seed 20261116

The input Parquet requires ``subset_id``, ``corrected_junction_nt_key`` and
``molecule_groups``. The processing TSV supplies ``input`` and ``molecule_groups``
for each library. A new output directory is mandatory. Compact curve/summary TSVs
and a source/parameter receipt are retained; replicate matrices are discarded.

Protein-level isotype evidence
------------------------------

``vdjaqp.isotype_evidence.isotype_evidence(placements, candidates, scores=None)``
derives isotype evidence from existing placements; it adds no search, scoring or
candidate sequence. A placement is **J–C spanning** when its interval contains an
FR4 residue (``region_ids`` 7) immediately followed by a constant residue (8). Only
constant positions a candidate already carries can be spanned: natively
alignment-supported contig positions (``contigs.adapt_contigs``) or residues of a
supplied reference protein. Reference C stubs are never appended, so variable-only
candidates yield no J–C evidence by construction.

The covered constant residues (the peptide's when a ``peptide`` column is given)
are compared under I/L equivalence with human IGHC CH1 references
(``IGHC_CH1``; source in ``IGHC_CH1_SOURCE``), indexed from CH1 residue 1. Per
identity class (``identity_class`` or I/L-canonical ``cdr3_seq``) the table
reports J–C placements and peptides, J–C peptides whose placements (in any
region) fall on that class only, reference-consistent peptides, the protein-level
isotype, a subclass only when every supporting peptide fits one gene, the
**identity-linked isotype** from class-unique J–C peptides alone (a shared
framework–constant peptide shows the isotype is present but does not tie it to one
V(D)J identity), the ``c_name`` isotype and agreement, and
whether any record passes the identity rule (``confirmed_nats >= threshold``,
default 20 nats). A CH1 prefix shared by several subclasses (e.g. IGHG1–4 over
residues 1–13) gives an isotype-level call. ``unresolved`` means consistent
peptides fit several isotypes; ``conflicting`` means several isotypes are each
specifically supported. J–C evidence links identity and isotype in one observed
peptide; it is not an abundance, isotype-fraction or subclass-titer estimate, and
heavy-chain references only are provided.
