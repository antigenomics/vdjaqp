Input and output contracts
==========================

Sequence-preserving AIRR-conditioned MS inputs
==============================================

``vdjaqp.airr_ms_inputs.prepare_airr_ms_inputs(candidates, *, physical_clone_ids, equal_protein_amounts, source_format, output, native_protein_intensity=10000.)`` prepares a fresh **private** bundle from an already annotated protein dataframe. It reuses the lossless adapter validation; coding nucleotides are optional. It does not generate sequences, reannotate, trim constant stubs, fill missing regions or infer physical membership.

The caller must supply a nonempty unique list of physical clone IDs and explicitly accept ``equal_protein_amounts=True``. These amounts are relative protein weights of 1.0 per selected record, not measured abundances or RNA counts. Invalid/noncanonical protein sequences and unknown/duplicate selections fail before output creation. Identical sequences with distinct candidate IDs remain distinct records. Existing output directories are rejected.

Returned paths:

- ``candidates``: unchanged original dataframe in compressed Parquet.
- ``search_fasta``: all input rows, with stable row-ordinal neutral IDs.
- ``proteins_fasta``: explicitly selected rows, in input order, with native MSSimulator intensity headers (default ``[# intensity=10000 #]``).
- ``id_mapping``: private original-ID/neutral-ID mapping for every row, physical membership and relative input amounts and native input intensities (null for absent physical proteins).
- ``manifest``: counts, provenance, explicit assumptions and artifact SHA256 hashes. The input hash identifies the unchanged dataframe snapshot, not the original upstream file serialization; callers retain original source hashes separately.

Both FASTAs preserve every supplied sequence byte. The mapping and candidate snapshot preserve original identities privately; neutral identifiers do not authorize sequence or artifact publication. These prepared inputs are requested artifacts of a fresh build, not reusable computed analysis caches. This helper runs no MS simulation, scoring or scientific analysis and chooses no acquisition scenario. The caller must obtain scientific population/amount/settings decisions before using it on experimental inputs.

The per-protein native intensity defaults to **10,000**, matching the `pinned MSSimulator input default <https://github.com/OpenMS/OpenMS/blob/60ce8966d577b43a6419403dec171539e08fe513/src/utils/MSSimulator.cpp#L203-L230>`__. Native intensity is not normalized away: digestion and ionization depend on its numerical scale. ``native_protein_intensity`` must be positive and finite, as must its population total. The manifest records both per-row intensity and total intensity in uncalibrated native input units, separately from equal relative weights of 1.0. Different physical population sizes consequently have different total native input intensities; this is not a measured experimental concentration or an abundance fit. Changing the default for a scientific run requires an explicitly justified setting.

Approved direct-MS driver
-------------------------

``scripts/simulate_airr_ms.py --repertoire /absolute/private/repertoire.parquet --target-ids /absolute/private/panel.csv --source-format 'actual supplied export provenance' --output /absolute/private/fresh-case --physical-scope panel --ms2-per-bin 5 --cpus 16`` runs one fresh case. The panel CSV must contain the twelve distinct supplied ``clone_id`` values. ``--physical-scope all`` instead includes every supplied protein physically; the search database always contains all candidates. The approved inventory is both scopes at capacities five and ten, fixed 1,200-second acquisition and two separate proteases per case. Outputs are private.

The driver reuses native simulation/search, spectrum-contributor linkage, the database-search error-model route and the unchanged scorer. Whole supplied/stub/raw-region/gene coverage is distinct from corrected panel-variable coverage. Entire peptides crossing the variable/stub boundary are excluded from the panel-variable comparison. Each native query batch includes all eligible called species and full candidate-variable evaluation queries; evaluation truth never changes inference inputs. Global identity classes are resolved before panel stratification. Spectra are removed only after successful validation/export; cleanup status records observed deletions.

``scripts/compare_airr_ms.py --targets /absolute/private/recorded-targets.tsv --case /absolute/private/completed-case --output /absolute/private/fresh-comparison`` accepts repeated ``--case`` flags. It reads final deliverables and the recorded targets, without rerunning simulation or historical scoring. It exports signed/absolute residuals for compatible panel-variable species and nine coverage regions, retaining exact denominators and shared-acquisition labels. It selects no weighted objective, tolerance, prior or posterior; recovery is excluded. Simulated total species include the constant stub and are not a matched experimental total-yield quantity.

After a case completes, ``scripts/compress_airr_ms_artifacts.py /absolute/private/completed-case`` losslessly packages large listed TSV, XML, JSON and FASTA artifacts. It validates report/cleanup status and hashes, verifies the gzip round trip before removing a plain copy, and records original-to-compressed hashes in ``retention.json``. Parquet, small summaries and the original result report remain unchanged. This is final-artifact storage, without scientific recomputation or intermediate reuse; repeated packaging is refused.

Fresh AIRR repertoire capture exports
=====================================

``vdjaqp.airr_capture.capture_tables(evidence, scores, reference_ids, ...)`` uses the existing encoded alphabet and ``expand_alignments`` in bounded chunks. It unions **agreeing** residues, including I/L equivalence, not every position of an approximate placement. Duplicate offsets/peptides never double count a residue. No per-clone native query, new search engine, posterior calculation or dependency is introduced. Junction classes are encoded once with class offsets; np.logical_or.at scatters agreeing support into class-position unions, and Polars group-bys form membership and count sums. Alignment expansion and union remain vectorized NumPy, the first sufficient optimisation rung.

Three views are explicit: all returned compatible alignments; alignments with existing llr>0; and all compatible alignments restricted to candidates at the existing confirmed_nats≥20 threshold. The last is not a posterior support cutoff. Neither observation ambiguity nor gene association proves unique source identity.

Six flat Parquet exports per context:

- ``records``: every record × view × supplied region, protein/junction interval and operational V/J/C region mask, retaining zeros and undefined fractions. Junction slices must exactly equal supplied cdr3_seq; anchored junction remains distinct from IMGT anchor-excluded CDR3. Zero C-length means unavailable supplied C region, not a zero-length physical constant.
- ``record_capture``: all/background/reference counts and fractions, residue totals, and optional background count-weighted any-capture/mean interval coverage.
- ``junctions``: exact anchored AA-string classes, distinct from record IDs and full variable mass classes; union coverage across member records, any partial and complete support, class membership counts. It does not merge I/L junction strings into a biological clone or claim unique recovery from partial support.
- ``mapping_ambiguity``: distinct compatible candidate-record and exact junction-string counts per peptide/view. Shared support remains explicit; no arbitrary winner.
- ``genes``: gene-call-bin totals and associated-record/direct-supplied-region capture. Allele suffixes use existing ``to_gene_level``. Compound calls remain native bins; they are not arbitrarily split into uniquely assigned biological genes. Empty, null and literal missing-call bins remain in the output and record denominators.
- ``gene_capture``: nonmissing called-bin presence/capture counts and fractions, excluding missing bins from the distinct-called-bin denominator; an empty called set has a null fraction, not a fabricated zero denominator.

The supplied dictionary is 0FR1,1CDR1,2FR2,3CDR2,4FR3,5CDR3,6NDN,7FR4,8C. ``supplied_V_framework_proxy_0_to_4`` uses labels0–4 only. Label5 includes both junction flanks and is explicitly excluded from this conservative V-framework proxy. ``supplied_J_FR4_proxy_7`` and ``supplied_C_stub_proxy_8`` are partial supplied proxies, not entire templated J or physical C sequences. Gene-associated capture (any supported record carrying the call) remains the primary called-gene metric. No full IMGT V/J/C placement is inferred. Label6 includes templated D, not pure N. Record summaries retain captured/all-record and captured/defined-record fractions separately, so missing C intervals are never silently treated as defined zeros.

Native ARDA protection masks
----------------------------

When ``junction_core_status`` is non-null, region 6 is a **computational junction protection mask**, not measured N additions or a biologically established NDN interval. Capture reports those rows as ``region=junction_protection``; inferred V/J-excluded cores can contain templated D, while unresolved geometry conservatively protects the entire anchor-excluded junction interior. Per-record exports retain ``junction_core_status``, ``junction_core_start_nt``, ``junction_core_end_nt`` and ``junction_core_protection_provenance``; bounds are inferred zero-based half-open junction-nucleotide offsets and remain null for unresolved geometry. Unknown masks must not be pooled or narrated as measured N coverage.

For these native rows, ``region=CDR3`` uses the explicit amino-acid interval ``[cdr3_start+1,cdr3_end-1)``, so protecting a residue does not remove it from the CDR3 coverage denominator. The anchored-junction view and exact-junction classes already use explicit coordinates and remain unchanged. Rows with null/absent native status, including supplied reference queries, retain their existing supplied-label masks and names. All-legacy calls keep their output columns and values unchanged; supplied label 6 retains its operational upstream meaning.

Optional counts are a DataFrame/TSV with exactly ``clone_id,count``, unique exact IDs covering **all background records and no inserted references**, nonmissing exact integer dtype, nonnegative counts and positive total. Floating storage is rejected, even for apparently integral values. Total counts must fit Int64 (≤2^63−1); negative/overflow/unsafe inputs fail. Nullable UInt64 per-record counts and exact integer aggregate sums preserve values beyond2^53. ``count_unit`` must be ``reads`` or ``UMIs``, with nonempty provenance; these units are not interchangeable. Missing counts fail. Inserted reference counts remain null, never zero. With supplied counts, weighted fractions are available for background only; mixed all-record and reference-only weighted totals remain null because their injected-reference amounts are unknown. Without supplied counts, no read-weighted values are invented. Experimental MS n_reads, synthetic RNA allocations and protein weights cannot supply this input.

Review-gated full export
------------------------

.. code:: sh

   .venv/bin/python scripts/profile_build.py --output private_data/airr-capture-resource.json --project-root . --workload-root private_data/airr-capture --timeout 1800 --rss-limit 24000000000 --disk-limit 10000000000 --cpus 16 -- .venv/bin/python scripts/export_airr_capture.py --input-root private_data/igseq_datasets --output private_data/airr-capture --cpus 16

Default aggregation is accepted corrected ``alignment_union``; ``--aggregation exclusive_spectrum`` is explicit for future authorised comparisons. The driver freshly calls Evidence.build for gel/solution and each original eight individually split target-excluded calibration contexts. No saved alignments/background/features are read. All twelve references define background/reference reporting populations, including individual runs where only one reference is the evaluated target. All source/raw input hashes, aggregation/count assumptions and output hashes are recorded; input/source drift fails. Fresh output is required. Private identifiers, junction sequences and raw mappings must never be committed publicly.

Engineering size: at33,052 records the record table has1,388,184 rows/context (three views × fourteen intervals), or13,881,840 rows for ten contexts. Compressed outputs are estimated at tens to a few hundred MB, not measured here; monitor actual disk/RSS before broader variants. The prior twenty-build correction took 524s; this is context, not a guarantee for ten builds plus capture exports. Starting qualification retains16CPUs/24GB/64threads/10GB and1800s guards. Two constructed checks cover shared partial junctions, agreeing-union duplicate/chunk invariance, positive/identified views and strict background read counts with unmeasured references. No actual cohort export was launched by the implementation worker.

Reviewed revision qualification: three constructed checks passed in0.18s. A representative constructed33,052-record,140-AA/record fixture with264,416 alignments completed capture in0.760429s (one Polars thread), with1,312,243,712B whole-process peak RSS on macOS. It produced1,388,184 record rows, 126 capture summaries and198,312 junction rows. Fixture construction is included in process peak memory but excluded from capture elapsed time. These are software qualification measurements, not real-cohort performance or CPU-scaling claims. An empty reference population retains zero-record/null-fraction summaries and no junction-class rows. Counts remain unavailable until an actual source is supplied; no request or invented surrogate count is added by this module.

All-locus observed compatibility capture
========================================

``scripts/export_all_locus_capture.py`` adds a descriptive capture report for the native IGH/IGK/IGL background. It does not alter canonical observations, native annotation, MS acquisition/search, scoring, references or identity criteria. Heavy scoring/recovery outputs remain separate. No light score or light identity is inferred.

.. code:: sh

   .venv/bin/python scripts/export_all_locus_capture.py \
     --geometry amplicon --observed-receipt OBSERVED/receipt.json \
     --native-clones NATIVE/smoke.clones.tsv \
     --acquisition-manifest SOLUTION_TRYPSIN/manifest.json \
     --acquisition-manifest HEAVY_BAND_TRYPSIN/manifest.json \
     --acquisition-manifest SHARED_CHYMOTRYPSIN/manifest.json \
     --output FRESH_REPORT --cpus 16

The same interface accepts ``--geometry rnaseq``. Exactly three completed same-geometry acquisition manifests are required, each pinned to the supplied canonical native ARDA observation receipt. The original solution trypsin, heavy-band trypsin and one shared chymotrypsin are read once each through the existing native idXML reader. The existing solution and heavy-band pools each include their own trypsin plus the shared chymotrypsin once. Pools share observations and are not independent replicates. A light-compatible peptide placement in a heavy-only acquisition does not imply an assayed light protein.

Population and source verification
----------------------------------

The universe is every valid native corrected-clone row retained in the canonical catalog, with explicit true background membership. Inserted reference rows are rejected, not silently incorporated. Canonical status, source format and geometry are required; no particular commit field or private-wrapper provenance is required. Catalog, heavy candidates and all-locus search entries are schema/hash/count/key/role verified. Optional producer/export commit metadata and producing code hashes remain provenance, without a fabricated version label. The supplied clone TSV must match the canonical producing ``.clones.tsv`` SHA256; native row ordinals join V/J/C calls and exact integer duplicate counts. Joined locus, junction, parent/source keys and counts must agree. Native source inventory reports retained/excluded rows and support by locus. Other native raw-TSV pins remain historical producer provenance and are not claimed reverified.

Representative native sequences come from existing search entries. Catalog rows legitimately absent from search entries are retained only when their source is a junction-only sequence shorter than the existing six-AA search minimum; their corrected junction remains the defined short target and receives zero compatible peptide capture. Longer missing entries fail. Sequences are not filled from truth or germline, and native normalization/correction provenance remains explicit. ``OBSERVED_SEGMENT`` means the available exported native segment, including explicitly labelled corrected-junction-only fallbacks; it does not imply complete variable sequence or uncorrected original bases.

One search, explicit coordinate unions
--------------------------------------

All distinct native target AA sequences are indexed once; aliases preserve every native record. All distinct called peptide AA sequences from the three acquisitions form one query batch through existing ``exact_peptide_placements``. Only the established six-AA minimum is applied; no PSM score/quality/identification threshold is added. Existing native idXML calls are compatibility evidence, not a new calibrated identification claim. Short queries remain in an exclusion table.

Condition selection and alias expansion occur through keyed native tables after the single search. Exact matching uses the established I/L equivalence. Native hits, alias/condition expansion and residue expansion retain explicit bounds: defaults5,000,000hit/association rows and2,000,000expanded slots per whole-chain chunk, as in existing evaluation helpers. Exceeding a bound fails; it does not truncate results. Native returned-hit and alias bounds are checked before subsequent expansion; the native matching API itself returns its hit table before that check.

``bounded_coverage``/``region_coverage`` compute residue unions separately for each native record and condition: available segment, conserved-anchor-inclusive JUNCTION and anchor-excluded IMGT CDR3. No region6mask is interpreted as N-only sequence. Exact anchored-AA junction classes are grouped within locus, with I/L-equivalent biological strings kept distinct. Clipped junction placements are unioned across class members; multiple partial peptides can jointly cover a complete class. Class support does not identify a unique biological clonotype. Mapping multiplicity records compatible native records, exact junction classes and loci; no shared peptide is assigned exclusively to one source.

Fractions and weighting
-----------------------

Record summaries give captured, partial-only and complete records over the declared native-record universe, with defined residue totals and coverage. Gene tables report **associated called-bin capture**, meaning any compatible peptide overlap with a record carrying that native call. They do not report direct V/J/C segment residue coverage. A scalar call has status ``single_call``; allele suffixes collapse to a ``call_bin``; raw calls remain in the record table. Compound/unresolved calls remain explicit bins without splitting or choosing a gene. Missing calls are separately counted. Native C calls can be isotype classes such as IGHG rather than genes or subclasses. Their bins carry ``native_C_call_not_gene_or_subclass_verified``; no C gene/subclass count is claimed. Raw native C strings remain unchanged. Metrics use ``distinct_single_call_bins``, not single-gene counts. Summaries expose both all called-bin fractions and the subset of unambiguous single-call-bin fractions.

The weight unit is ``native_processed_fragment_support``, supplied by native corrected-clone ``duplicate_count`` after fastp merging/correction. For a supported record, its count enters a numerator once, regardless of placement/alias multiplicity. Record-weighted capture divides supported native counts by total counts in the same native-record population. Weighted mean coverage divides summed count×coverage by native counts for defined intervals. Zero denominators remain null. These counts are neither original raw pairs, measured RNA abundance nor experimental AIRR reads.

Junction-class summaries aggregate member counts once per exact-AA class. Their weighted capture describes **member support of a class with union capture**; it does not mean every member independently had a compatible placement. It remains separate from record-weighted coverage. ``input_totals.parquet`` separately records original raw input pairs, fastp processed-fragment records and retained catalog native support; none is substituted for another denominator. Fastp pair/merged/unmerged arithmetic must agree.

Observation-only repertoire export
==================================

Use the completed assembly directory as ``--root``; this command performs no read simulation, assembly, annotation regeneration, truth lookup or reference insertion:

.. code:: sh

   .venv/bin/python scripts/export_observed_repertoire.py \
     --root /path/to/completed/run --geometry amplicon \
     --output /path/to/fresh/export --cpus all

Inputs are ``assembly/raw_support.json``, ``assembly/arda_contigs.json`` and ``assembly/read_support/{geometry}_support.json``. Raw/native rows use ``sequence_id`` and ``sequence``. The selected geometry’s support rows use ``contig_id``; their keys become ``{geometry}_{contig_id}``. The native adapter validates sequence equality and retained-primary exact-match masks. Geometry selection uses the assembly’s existing mode prefix, never a truth identifier.

For transcript-aware RNA, pass ``--geometry rnaseq --boundary-contexts`` with an explicit public/model context JSON array containing ``target``, ``sequence`` and ``mature_start`` (1-based closed nucleotide coordinates). The existing boundary search runs one internally threaded batch. Without this option the native adapter runs without mature-boundary projection, as in the original coding-only branch; this is recorded explicitly. No generating-chain offsets are accepted. The raw preimputation consensus is retained unchanged even when boundaries are projected.

The fresh output contains ``catalog.parquet``, ``candidates.parquet``, and ``search_entries.parquet``, plus selected raw/native/support JSON, optional boundary and alignment outputs, and ``receipt.json``. The receipt records exact input/code SHA-256 hashes, source paths, resolved CPU budget and actual counts by observed locus. Native candidate construction remains heavy-only; light observations and supported segments remain in catalog/search entries. No fixed heavy count is imposed. Every exported observation has ``background_member=true`` from its source role. Query insertion belongs to ``vdjaqp.reconstruction.append_known_references``; this CLI neither inserts queries nor changes native candidate features. Simulation truth may be joined separately for evaluation only. Existing output directories are rejected, including incomplete earlier exports.

Evaluation-only scientific metric exports
=========================================

``vdjaqp.simulation_metrics`` reuses Polars vectorized residue expansion/joins and one batched native ``PeptideMatcher`` call. No new dependency, mapper, simulation model, threshold or scoring route is introduced. Truth inputs are for evaluation only; this module never feeds them into the inference model.

- ``exact_peptide_placements(chains, peptides, threads=None)`` takes ``chain_id, protein_aa`` and ``peptide_id, seq``. It returns ``chain_id, peptide_id, start, end`` with zero-based half-open protein coordinates, all repeated exact placements and I/L equivalence. Null/duplicate identifiers, noncanonical sequences and queries shorter than the existing six-residue search minimum fail explicitly rather than disappear from a coverage denominator. Empty typed inputs return typed empty placements. One existing native batch owns searching; its resource resolver owns clamping.
- ``region_coverage(regions, placements, evidence_levels)`` takes unique ``chain_id, region`` rows with ``locus, start, end`` and optional scalar gene metadata. Placement rows contain ``chain_id, start, end, evidence_level``. All region metadata survive. Input ``status`` is retained as ``region_status`` (for example ``not_applicable`` for light-chain D); output ``status`` describes defined/undefined coverage. Supplying both status names is rejected as ambiguous. Explicit evidence levels preserve zero-covered regions even when no placements exist. Results add integer ``region_length, covered_residues``, fractional ``coverage``, and ``status``. Undefined regions must be explicit rows with both endpoints null; zero-length regions also have undefined coverage. An absent locus cannot be inferred from missing rows. Full true coordinates supply denominators; partial candidate lengths must not replace them. The region set must include the chain’s full coordinate bound, otherwise placements beyond its declared maximum fail. Overlap/repeated hits are deduplicated by residue before summation; named overlapping views such as junction remain separate. No macro/pooled averaging is silently applied.
- ``psm_correctness(peptides, native_truth)`` takes the native peptide-export frame (``peptide_id, seq, native_spectrum_reference`` required) and the native MS2 contributor JSON list. It returns ``(per_psm, summary)``. Per-PSM rows retain input order/metadata and add nullable ``correct``, ``status``, integer ``native_contributor_count, matching_protein_count``, and nullable ``protein_ambiguous``. A match to any I/L-equivalent contributor is correct. Protein compatibility is not unique clone recovery. Missing truth links are explicitly unscorable; missing spectrum IDs, duplicate truth IDs and noncanonical modified sequences fail. Native contributor fractions are not assignment probabilities.

The flat summary has ``quantity, count, denominator, fraction, unit``. Count cells are exact integers; for the ``correct_psm_rate`` row, ``count`` is the correct-assignment numerator, ``denominator`` is evaluated assignments and ``fraction`` is their ratio (null at zero denominator). Native event count, called/evaluated spectra, correct-call spectra and PSM assignment counts are separate quantities. No unlinked assignment enters the accuracy denominator. No rank/recovery cutoff is introduced.

Five focused regressions cover repeated/I-L placements and explicit short-query failure; union/zero/undefined/out-of-range coverage; coisolated contributor sets, wrong/unlinked calls, empty exports and duplicate keys; and typed zero-hit coverage/missing spectrum IDs. They were observed failing before the corresponding implementation/fixes. Verification command:

A schema/compatibility inspection of the retained native software fixture accepted 15 of 15 PSM assignments in each scenario against 154 native MS2 contributor sets. This was evaluation of retained fixture exports, not a new simulation, real-data rerun, fitted accuracy estimate or scientific benchmark. Parent integration must produce region/truth tables and call these APIs during its fresh run before releasing intermediate alignments.

Heavy-variable identity coverage and recovery
---------------------------------------------

``identity_metrics(metadata, candidates, ranks, associations)`` returns ``(identity_population, candidate_identity, summary)``. Metadata requires physical ``chain_id, locus, variable_aa, protein_aa``; optional ``variable_length`` must agree with the full variable sequence. Candidate keys and ranked keys must match exactly; ranks require integer positive deterministic ``rank, n_tied`` and finite ``confirmed_nats``. Associations are exact evaluation-only fragment placements with ``clone_id, chain_id, start, end``; unknown keys, noninteger/invalid/out-of-protein bounds fail. Native scoring remains unchanged.

The inference endpoint is the **heavy-variable amino-acid mass-equivalence class**, defined by I→L normalization of the complete variable sequence. A deterministic opaque SHA256 class key groups physical chains sharing that sequence, without claiming unique nucleotide-clone, lineage, constant-isotype or paired-antibody identity. Every physical true-heavy row remains in the population table, including absent identities. Class-level summaries deduplicate those rows and multiple candidate aliases.

An association is full-variable compatible only when it starts at protein position zero and reaches at least ``variable_length``. Partial exact fragments establish compatibility, never full-variable identity. Within each candidate/class pair, a full association supersedes partial associations; repeated placements cannot inflate candidate counts. Per-candidate ``full_class_count, partial_class_count`` distinguish ``unique_full``, ``ambiguous_full``, ``partial_only``, and ``no_exact_compatibility``. Missing exact compatibility can reflect assembly errors or partial input and is not automatically a false identity.

The fixed accepted criterion is **confirmed_nats ≥ 20**, including its boundary; combined score is not substituted. No threshold scan or new framework-only filter is added. A mass class is recovered only if at least one shortlisted candidate fully matches exactly that one class. Full matches spanning several variable classes remain ambiguous. Population fields include full/partial-compatible candidate counts, best full-compatible candidate ID/rank/tie count (nullable), and ``class_recovered_at20``. Best rank retains the scorer’s deterministic order and tie count and does not by itself prove identity resolution.

Summary quantities distinguish physical true-heavy rows, true variable mass classes, full-variable-observable classes, recovered classes over **all** true classes, and conditional recovery over full-variable-observable classes. A zero conditional denominator gives null recovery. Shortlisted candidates, uniquely full-compatible shortlisted candidates and unresolved shortlisted candidates remain separate integer counts; no fabricated precision estimate is reported. ``count`` is the numerator in fractional summary rows. Input candidate frames should not contain rank/evidence columns; the separate ``ranks`` argument is authoritative.

Seven focused tests now pass, including same-class physical duplicates/I-L equivalence, repeated candidate aliases, full versus ambiguous partial compatibility, absent identities, exact 20-nat boundary, ambiguous full-variable matches, empty populations, missing rank rows and nonfinite evidence. Compatibility inspection accepts the retained fresh-five table schemas without a native rerun. Full scientific population/recovery results must come from the parent’s fresh integrated run.

Explicit reference insertion for matched reconstruction
=======================================================

The author-approved matched design adds 12 known reference rows to a FASTQ-derived background. This construction is separate from the observation-only generic simulator. ``vdjaqp.reconstruction.append_known_references(candidates, search_entries, reference_panel, expected_count=12, source_format=...)`` performs that explicit append; it neither generates nor annotates sequences. Both background tables must be unmodified ``adapt_contigs`` outputs. The supplied heavy panel must already contain scorer columns and coding ``nt_seq``; the existing AIRR-compatible adapter checks its sequence/junction/region schema. Additional validation requires canonical protein/coding DNA, frame-zero translation agreement, conserved junction anchors, heavy V/J calls, unique IDs, and no collisions with either background table.

The resulting ``candidates`` table preserves all original background fields and adds ``candidate_construction=fastq_derived``. Inserted rows use ``candidate_construction=reference_panel_insertion``, ``sequence_origin=supplied_known_reference``, and ``coverage_kind=supplied_reference_sequence_not_rna_support``. Their contig support/coordinates are null, not fabricated. Panel sequences also enter the explicit matched search database with their own origin label. The returned ``reference_panel`` is the supplied panel in AIRR-compatible interchange format, preserving its actual producer provenance and extra source columns. Input identifiers are retained. Identical sequences under different supplied candidate IDs are not merged; downstream mass-equivalence evaluation owns that separate denominator.

The function returns three frames and modifies neither input. The generic ``smoke_metrics.py`` driver does not call it. No evaluation-truth object, automatic planted-reference selection, germline completion, reverse translation or scorer change is introduced.

Standalone construction
-----------------------

.. code:: sh

   .venv/bin/python -m vdjaqp.reconstruction \
     --candidates /absolute/background_candidates.parquet \
     --search-entries /absolute/background_search_entries.parquet \
     --reference-panel /absolute/supplied_known_panel.parquet \
     --expected-count 12 --source-format 'actual known-panel annotation producer' \
     --output /absolute/fresh/matched-inputs

The absent output directory receives three Zstandard Parquet files and ``insertion_manifest.json``, with role-keyed input/output SHA256 hashes and exact row counts. Input validation finishes before directory creation. Existing output is rejected, including a competing directory creation, rather than overwritten. No files are uploaded. This is a construction step in a fresh matched build, not authorization to reuse a prior simulated background or resume completed stages.

Four focused tests pass in 0.30 s: append and preserve background at panel sizes one and twelve; reject wrong counts, collisions, invalid coding/heavy calls and repeated insertion; and CLI count/hash/no-overwrite behavior. The twelve-row fixture deliberately repeats one sequence under distinct IDs to verify candidate-row retention; it is not a simulated biological panel. No experimental data, scoring, assembly or scientific parameter run is used by these tests.

Remaining integration: generate the full synthetic matched population and a separately annotated known synthetic panel, apply the declared supplied constant-stub policy, then invoke this step in the fresh matched workflow. The physical protein pool remains separate from the RNA background. This adapter does not establish an observed background count of 33,040, model real isotype/SHM distributions, fit assay settings or demonstrate biological recovery.

Evaluation-only chain and residue-region exports
================================================

``vdjaqp.simulation_truth.annotate_truth(chains, threads=None)`` calls installed ARDA ``reannotate_contigs`` once on every complete generated variable nucleotide sequence and returns ``(chain_metadata, regions)`` as Polars frames. Empty input returns typed empty frames. The CPU budget uses the existing shared resolver; no per-chain calls, pool or new dependencies were introduced. Existing batched native annotation suffices.

For a run retaining full native annotations, make that same one NT batch in the driver, persist the returned native records, and call ``regions_from_annotations(chains, annotations)``. Do not call both routes and annotate twice. Neither output belongs in scorer inputs.

Metadata preserves generator V/D/J calls, physical C gene/accession, family/parent/antibody keys, RNA weights, secretion factors, mutation-event and observable-substitution counts separately. Annotator locus/V/J calls occupy separate fields and never replace the generating model’s calls. Full physical protein and variable sequences/lengths retain their original denominator. Supplied complete coding RNA is checked against physical protein plus terminal stop; a native untranslated terminal partial codon is retained explicitly, not silently used as an amino acid or removed from provenance.

Each region row has ``chain_id``, ``locus``, ``region``, ``start``, ``end``, ``status``, ``source``. Coordinates are zero-based half-open **amino-acid** intervals. ARDA annotation NT endpoints are one-based inclusive; GenAIRR genotype endpoints are zero-based half-open. Out-of-range/reversed endpoints, invalid key joins, reverse coding orientation, overlapping partitions or sequence disagreement fail explicitly. FR4 ending at the native terminal partial codon projects only its complete variable codons. Non-codon-aligned boundaries otherwise retain null endpoints and ``unavailable_non_codon_boundary``. FR1/CDR1/FR2/CDR2/FR3/CDR3/FR4 are annotation-derived, with the full physical constant suffix ``C`` appended separately. A complete partition must cover the whole variable sequence exactly. Missing native boundaries produce null endpoints and ``unavailable``, never fabricated positions. No NDN partition is inferred; the separate D view described below uses native genotype segment coordinates, not a non-templated designation.

``JUNCTION`` is an overlapping known-truth view including conserved Cys and terminal Phe/Trp. ``GENERATOR_CDR3`` is the separate overlapping known-truth view excluding both anchors: ``[junction_start_aa + 1, junction_end_aa - 1)``. Both use source ``generator_junction``. Neither replaces annotation-derived ``CDR3``, and neither may be summed into the disjoint FR/CDR partition. Use ``GENERATOR_CDR3`` when the desired denominator is the known generated anchor-excluded CDR3; use ``CDR3`` to inspect the inferred annotation and its availability.

Native genotype fields ``junction``, ``junction_aa``, ``junction_start`` and ``junction_end`` are required. Starts/ends must be integers, zero-based half-open NT coordinates on complete codons. The physical NT slice must exactly equal native ``junction`` after case normalization; its translation must equal ``junction_aa``, and its AA position must equal the uniquely placed anchored junction in the physical variable protein. Invalid slices, translations, frames or coordinates still fail. No coordinates are manufactured from ARDA for older records missing these fields.

Metadata retains the original native junction string, AA sequence and NT endpoints as ``generator_junction_nt``, ``generator_junction_aa``, ``generator_junction_start_nt``, and ``generator_junction_end_nt``. Raw annotation junction strings occupy ``annotated_junction_nt`` and ``annotated_junction_aa``. ``junction_annotation_concordance`` is typed string metadata: ``matches``, ``missing`` (both inferred strings absent), or ``conflict`` (any supplied inferred string differs). An alternative inferred J alignment does not invalidate sequence-verified known truth or replace generating V/D/J calls.

FR3/CDR3/FR4 inferred boundaries are checked against known anchors: Cys belongs to FR3 and terminal Phe/Trp to FR4. Conflicting boundaries stay null with ``unavailable_junction_boundary_conflict``; coordinates are never shifted or filled from germline. These nulls are applied before checking whether an available complete partition covers the variable sequence. Out-of-range/reversed boundaries and overlaps still fail. Every physical chain is retained, including chains with missing/conflicting annotation; known ``JUNCTION``, ``GENERATOR_CDR3`` and physical constant denominators remain available independently of those inferred masks.

API correctness evidence
------------------------

On 2026-11-10, one diagnostic AA batch and one NT batch were run on the same already retained six-chain synthetic fixture, each using 16 threads; no generation, read simulation or historical analysis was rerun. The AA call returned six rows in 1.006 seconds but projected a heavy native junction one amino acid late (true C…W became A…G). It is therefore not used. This observation diagnoses an API-boundary mismatch, not its upstream source-code cause; no ARDA source was changed.

The NT call returned six rows in 1.009 seconds. All six native junction amino-acid strings matched generator truth exactly; conversion produced six metadata rows and 54 region rows (nine regions per chain, including overlapping junction), all available. Diagnostic raw native annotations are ``/tmp/truth-nt-annotation.json``; fresh scientific/integration runs must generate and retain their own native annotations rather than reuse this diagnostic.

Five focused tests pass, following an observed missing-module failure before implementation: three locus-schema/partition/anchor/constant cases, missing or invalid native bounds and protein boundaries, and typed empty/full-RNA validation. Actual NT diagnostic rows also pass the same validator. This verifies schema/coordinate handling on the tiny fixture, not population realism or comprehensive lambda execution.

Native genotype gene-segment views
----------------------------------

Region output additionally includes overlapping ``V``, ``D`` and ``J`` segment views and a scalar ``gene_call`` column (``C`` uses ``constant_nt_gene``; FR/CDR/junction rows have a null gene call). Gene and allele aliases, including ambiguous J2/3 labels, remain unchanged. These views are not an additional disjoint partition and must never be summed with FR/CDR or junction coverage.

The inspected GenAIRR source assigns segment starts/ends from native ``range.start.index()``/``range.end.index()`` in ``engine_rs/src/airr_record/builder.rs``; Python outcome export passes those integers unchanged. Its ``align_to_germline`` API explicitly slices ``sequence[v_sequence_start:v_sequence_end]``. These are zero-based half-open NT intervals, despite the AIRR-style field names; they are not converted using the ARDA one-based convention. Only fully templated codons are retained: AA start is ``ceil(nt_start/3)`` and AA end is ``floor(nt_end/3)``. Mixed-boundary codons are excluded, not counted as V/D/J sequence. Native D remains templated, never relabelled NDN. Missing segments or segments containing no whole codon have null endpoints; light-chain D is explicitly ``not_applicable``, not measured zero coverage.

The extension adds a sixth focused regression covering partial codons, gene-alias preservation, physical C and missing light-chain D. Revalidating the existing six NT diagnostic rows without another annotation call produces 72 region rows (12 views per chain). Heavy D is available in three rows; light D is not applicable in three rows. These annotation and segment populations have distinct denominators.

Native mutation counters remain unchanged source fields, not a validated substitute for cumulative substitutions from a known ancestor. In parent-state continuation, ``n_observable_mutations`` must not be assumed to equal the descendant-to-root difference count. Scientific SHM-family summaries still require explicit ancestor/event semantics; the current fixture verifies inherited sequence changes and retains the raw counters without relabelling them as germline mutation frequencies.

The 240-chain resource fixture exercises partial/incomplete annotations at population scale; its exact missingness counts and whole-build verification are in ``methods guide``.

Exact IGH V reference matches now appear as additional scalar chain metadata; native gene calls remain unchanged. Unique, ambiguous, unmatched and unevaluated states are distinguished without expanding chain populations. See ``methods guide``.

The truth-adapter regression now includes a software-string junction containing an internal W with an earlier alternative inferred J call. It preserves the known full junction and anchor-excluded CDR3, retains the raw conflicting annotation, and nulls only the conflicting inferred CDR3/FR4 intervals. Additional checks reject corrupted native junction NT and shifted or missing native coordinates. This is adapter verification without sequence generation, reannotation or scientific reruns. Historical fixture row counts above describe their original export schema; current exports add ``GENERATOR_CDR3`` as an additional overlapping view per chain.

A codon-aligned annotation can still cover only part of the physical variable: for example, a mutated leading codon may leave the inferred V/FR1 alignment starting at nucleotide 4. Valid FR1 intervals starting after AA zero and valid FR4 intervals ending before the physical variable end are retained as unavailable outer-region annotations: endpoints become null with ``unavailable_incomplete_variable_partition``. No coordinate is shifted and no sequence is filled from germline. Physical variable/protein sequence and lengths, known ``JUNCTION``/``GENERATOR_CDR3``, and the physical constant boundary remain intact. Out-of-range/reversed intervals, overlaps and internal gaps in an otherwise complete inferred partition still fail; a partial outer interval does not excuse an internal gap. Native annotation records retain their originally inferred bounds.

Explicit unit-count placeholders
================================

The author-requested interim publication recount uses ``count=1`` for every
candidate, including inserted references. ``capture_tables`` requires
``count_unit='placeholder'`` and explicit provenance for this mode, validates
complete candidate membership and exact integer ones, and adds
``counts_are_placeholder``/``count_unit``/``count_provenance`` to its outputs.
Measured reads/UMIs retain the strict background-only contract and leave
inserted-reference counts unmeasured. Unit-weight totals equal record totals;
constant counts cannot support a count–peptide-count correlation.

``vdjaqp.report.peptide_support_table(evidence)`` counts distinct input peptide
rows, I/L-equivalent species and zero-substitution species from the fresh native
match table, preserving zero-hit candidates. Shared species remain compatible
with multiple candidates and are not treated as unique identities.

Run ``./analyze.sh`` (or ``make analyze``) for the complete ten-context
experimental recount. Its source/raw/output hashes and process-tree resource
report accompany the completed private output. Failed runs preserve the previous
completed result. This is not a simulation-model or abundance-fit rerun.

Native junction eligibility
---------------------------

Observed-frame projection preserves ARDA's existing complete-junction flag.
Only complete, in-frame anchored junctions can receive a V/J-excluded protection
mask. An incomplete observed junction remains a partial observation; the exporter
never invents its conserved anchors. Native parent export prefers a complete
observed frame and otherwise declares a corrected-junction-only source, retaining
counts of complete and partial matching observations. Compatibility reporting can
retain partial protein intervals separately from complete-clonotype identity
queries. This is an annotation-domain check, not a read-depth or singleton filter.
