Data and provenance
===================

The companion data are maintained in the Hugging Face dataset
``isalgo/igprotseq``. Access follows the dataset settings; run ``hf auth login``
when authentication is required. Every file is
listed in its ``MANIFEST.tsv`` with size, SHA-256, group and local source path.

Bootstrap
---------

.. code-block:: sh

   uv sync --extra data
   .venv/bin/vdjaqp-data                                # calibration only (default)
   .venv/bin/vdjaqp-data --group donor simulations      # one-donor AIRR/peptides, simulations
   .venv/bin/vdjaqp-data --group external               # processed public gut cohorts
   .venv/bin/vdjaqp-data --group results --dest ../2026-vdjaqp-ms   # paper result tables
   .venv/bin/vdjaqp-data --verify DIR --group all       # re-check an existing copy

Each download is verified against the manifest SHA-256. Pin a dataset commit with
``--revision`` for exact reproduction. Mutable revisions are resolved once to an
immutable dataset commit before fetching the manifest and files. Manifest paths,
required metadata, sizes and digests are validated before publication.

Writers serialize through a bounded destination lock. Destination file identity
is checked to reject aliases of that lock before publishing any asset. Each file is copied to a
same-filesystem staging file, verified and atomically replaced; the manifest is
published last. This is per-file publication, not an atomic multi-file snapshot
for readers. Finish bootstrap before starting consumers. Concurrent local writers,
interrupted copies and killed lock owners are tested; network-filesystem locking
has not been qualified. Abrupt termination may leave an unused staging file.
``--verify`` and ``--dry-run --manifest PATH`` operate without network access.

.. list-table:: Dataset groups
   :header-rows: 1

   * - Group
     - Contents
     - Size
   * - ``calibration``
     - frozen error model (confidence knots [0, 1] → error [0.10, 0.01], ``mass_prior``
       substitution), robustness spec, declared method/design/wet-lab/OAS parameters,
       ARDA reference pin
     - < 1 MB
   * - ``donor``
     - access-controlled annotated repertoire and peptide example inputs;
       access and redistribution permissions are recorded with the dataset
     - 0.4 GB
   * - ``simulations``
     - three joint-workflow populations; simulation result tables
     - 0.1 GB
   * - ``external``
     - Crohn GSE261086 (RNA-seq AIRR, expression) and PXD050335 (tissue MS proteins,
       peptides, placements, VDJ-AQP scoring); coeliac MSV000081328 and PXD013777;
       PROTECT GSE109142 AIRR and outcomes; gene matrices
     - 2.8 GB
   * - ``results``
     - bulky manuscript result tables, same relative paths as the paper repository
     - 0.5 GB

Local data storage
------------------

``VDJAQP_DATA`` selects the local data root (default ``./data``, shared with
``vdjaqp-data``); :mod:`vdjaqp.paths` exposes
``DATA``, ``EXTERNAL`` and ``GUT``. Keep raw inputs read-only and retain processing
receipts with each result. Host-specific storage mappings and scheduler accounts
belong in private operational records, outside the distributable documentation.

Software versions
-----------------

``uv.lock`` pins the Python environment. ``vdjtools`` moved from 5.0.0 (removed from
PyPI) to 5.9.2 on 18 November 2026; IGH Pgen was bit-identical on 3,000 donor
junctions. Native tools on the cluster are pinned by SHA-256 in their installation
receipts.
