Parameters
==========

Defaults are read from ``nextflow.config``. Paths are strings; ``null`` means
the corresponding branch requires an explicit path. Numbers below are not fitted results.

.. list-table:: Workflow parameters
   :header-rows: 1
   :widths: 20 30 50

   * - Parameter
     - Default / type
     - Meaning or units
   * - ``mode``
     - ``'analysis'`` / string
     - branch selector
   * - ``samples``
     - ``null`` / path/string
     - analysis CSV
   * - ``input_root``
     - ``null`` / path/string
     - original input bundle directory
   * - ``antibodies``
     - ``1200`` / integer
     - paired antibodies
   * - ``seed``
     - ``20261110`` / integer
     - generation RNG seed
   * - ``read_seed``
     - ``2468`` / integer
     - native read-error RNG seed
   * - ``memory_fraction``
     - ``0.5`` / float
     - fraction of related SHM-family records
   * - ``mutation_rate``
     - ``0.01`` / float
     - mutation probability per eligible nucleotide
   * - ``protein_concentration``
     - ``'equal'`` / string
     - equal or head80 abundance allocation
   * - ``pairs_per_chain``
     - ``100`` / integer
     - read pairs per chain
   * - ``gradient_seconds``
     - ``1200`` / integer
     - seconds
   * - ``stage_timeout``
     - ``900`` / integer
     - seconds per stage
   * - ``arda_overlay``
     - ``"${projectDir}/.tool-envs/arda-overlay-read-map/overlay"`` / path/string
     - pinned native ARDA overlay directory
   * - ``repertoire``
     - ``null`` / path/string
     - as-is candidate Parquet
   * - ``panel``
     - ``null`` / path/string
     - CSV with exactly 12 distinct inserted-reference clone_id values
   * - ``references``
     - ``null`` / path/string
     - checksum-pinned public raw references
   * - ``additional_references``
     - ``null`` / path/string
     - additional public raw constant-reference directory
   * - ``producer``
     - ``'supplied annotated export'`` / string
     - actual source producer label
   * - ``physical_scope``
     - ``'all'`` / string
     - all or panel; does not alter search candidates
   * - ``ms2_per_bin``
     - ``5`` / integer
     - MS2 acquisitions per simulator bin
   * - ``outdir``
     - ``'results/workflow'`` / string
     - fresh published result root
   * - ``python``
     - ``"${projectDir}/.venv/bin/python"`` / path/string
     - project Python executable
   * - ``cpus``
     - ``Runtime.runtime.availableProcessors()`` / integer
     - CPU cores

Validation and interpretation
-----------------------------

``cpus``, ``antibodies``, ``pairs_per_chain`` and ``ms2_per_bin`` must be positive integers.
Seeds must be nonnegative integers. Numeric fractions, rates and times must be finite;
the producing Python APIs additionally validate biological ranges and family divisibility.
The population CLI restricts memory fractions to 0.25, 0.5 or 0.75 and read seeds
to integers below 2**32. Direct ``smoke_metrics.py`` defaults differ from the
workflow: its gradient is 120 seconds; Nextflow explicitly passes 1,200 seconds.
Supplied-AIRR branches accept MS2 capacity 5 or 10 only.
``memory_fraction`` describes composition, not measured cell phenotype or abundance.
``head80`` allocates 80 percent of each component protein mass to exactly
10 percent of its families, balanced within component/light strata; family counts
must be divisible by 20. RNA read pairs remain equal per observed chain.
It is a configured allocation, not an inferred experimental quantity.
Physical substrate membership, search membership and background membership remain independent.
CPU requests are capped to the actual allocation; CPU count enters native acquisition provenance.
