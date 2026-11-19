Tutorial: shared versus distinguishing peptide evidence
=======================================================

This five-second example uses only invented sequences, so it is safe to share and
needs no optional tools. It runs the same search, background, error-model and
aggregation path as the workflow (``Evidence.from_frames``).

.. code-block:: sh

   .venv/bin/python examples/toy/run_toy.py --output .build/toy-example

Inputs
------

``examples/toy/run_toy.py`` builds 303 heavy-chain candidates on one IGHV3-23
framework: a **target**, a **junction relative** differing from it by one junction
residue, a **CDR2 variant** with a different CDR2 and junction, and 300 background
sequences with random junctions. Five tryptic peptides come from the target: three
framework peptides shared by almost every candidate, one crossing into the junction
start and one spanning the junction interior into FR4.

Candidates need ``clone_id``, ``seq``, per-residue ``region_ids``
(0 FR1, 1 CDR1, 2 FR2, 3 CDR2, 4 FR3, 5 junction anchors, 6 junction interior,
7 FR4, 8 C), anchor-inclusive junction coordinates ``cdr3_start``/``cdr3_end`` and
V/J/C calls. Peptides need an identifier, sequence, per-residue confidence (``NaN``
for database-search calls), search type, de novo confidence and rank, and a read
count; see :doc:`input-output`.

Result
------

.. code-block:: text

   clone_id            peptide_matches  confirmed_nats  junction_nats  accepted
   target              5                26.3            25.7           true
   junction_relative   5                 2.9             2.8           false
   background_114      4                 0.13            0.13          false
   ...                 4                <0.02           <0.02          false

Every one of the 303 candidates matches at least one peptide, so a peptide-count rule
cannot separate them. VDJ-AQP credits each residue by its rarity in the candidate
background and lets candidates compete for each observation: shared framework
residues contribute almost nothing, while the junction peptide confirms 25.7 nats for
the target. The junction relative also matches all five peptides but receives under
3 nats, because the residue that separates it from the target contradicts the
observed peptide. Only the target passes the fixed 20-nat identification rule.

``framework_only`` flags candidates with at least one confirmed nat but none inside
the junction; such calls remain explicitly ambiguous. Change ``PEPTIDES`` in the
script (for example, drop the junction peptide) to see the target lose its
distinguishing evidence.

From example to cohort
----------------------

The same scorer serves the Nextflow workflow (:doc:`pipeline`): per-sample candidate
and peptide Parquet tables, an independently calibrated error model and fresh
output directories. Supported inputs include UMI-tagged amplicon and fragmented
RNA repertoires annotated by ARDA, database-search and de novo peptide calls,
separate or pooled protease acquisitions, and heavy, kappa and lambda chains.
Known-truth simulation branches generate repertoires, reads and spectra for
benchmarking new assay designs.
