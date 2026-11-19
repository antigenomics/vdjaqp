.. image:: _static/vdjaqp.png
   :alt: VDJ-AQP logo: skull with branching antibody antlers, peptides and mass-spectrometry rays and peaks
   :width: 320px
   :align: center
   :class: dark-light

VDJ-guided Antibody Quantification from Proteomes
=================================================

Reproducible antibody proteomics–repertoire attribution and simulation.

This package implements VDJ-AQP, the method of the accompanying manuscript; the Python package
and CLI are ``vdjaqp``. It links LC-MS/MS antibody peptides to B-cell receptor repertoire candidates.
Candidates compete for every peptide observation, residues are credited by their
rarity in the repertoire, and results keep compatible placement, credited evidence
and accepted identity separate. Start with :doc:`tutorial`.

* **Repertoire inputs:** UMI-tagged amplicon or fragmented RNA reads processed by
  assay-appropriate MIGEC molecule processing where supported and ARDA annotation
  for both sequencing geometries; subsets and their qualified pooled union in one candidate universe.
* **Peptide inputs:** database-search and de novo calls; trypsin, chymotrypsin and
  other proteases, separately or pooled.
* **Loci:** heavy, kappa and lambda chains, with region-, gene- and junction-level coverage.
* **Outputs:** per-candidate confirmed and junction evidence, framework-only flags,
  keyed placement tables and repertoire capture summaries.
* **Isotype evidence:** peptides spanning the FR4–constant boundary of candidates with
  native constant residues give protein-level isotype calls checked against CH1
  references and ``c_name`` (:doc:`methods`).
* **Simulation:** known-truth repertoires, reads and spectra (GenAIRR, InSilicoSeq,
  OpenMS) for benchmarking assay designs and comparators.

Choose the documentation for your task:

.. toctree::
   :maxdepth: 2
   :titlesonly:

   tutorials
   how-to
   reference
   explanation
