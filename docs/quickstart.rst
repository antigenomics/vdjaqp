Quickstart
==========

Run from the checkout root with Python >=3.11, Java 17–25 and Nextflow 25.10.2.

.. code-block:: sh

   sh setup.sh
   .venv/bin/python scripts/check_environment.py
   .venv/bin/python -m pytest tests -q
   nextflow run main.nf --samples samples.csv --outdir results/new-analysis

Use ``workflows/samples.example.csv`` as the neutral sample-sheet template. Inputs must exist;
output roots must be fresh. Scientific calculations require an approved input and
parameter specification. No resume, cache or generated intermediate reuse is supported.

Observed-analysis resource profiles
-----------------------------------

Use ``-profile analysis_node`` inside a scheduled node allocation, or
``-profile analysis_slurm`` for native Slurm submission. Declare the total
``--cpus`` and ``--analysis_memory_budget``; task defaults are one CPU and 4 GB.
Concurrency is bounded by both budgets. The Slurm controller CPU allocation does
not limit submitted task requests. These profiles support observed analysis only;
raw-data preprocessing uses dataset-qualified Slurm recipes.

::

   nextflow run main.nf -profile analysis_node --samples samples.csv --cpus 16 --analysis_memory_budget '64 GB' --outdir results/fresh-node
   nextflow run main.nf -profile analysis_slurm --samples samples.csv --cpus 16 --analysis_memory_budget '64 GB' --analysis_partition medium --outdir results/fresh-slurm

Optional native tools
---------------------

Install the committed recipes in ``envs/generation``, ``envs/preprocess`` and
``envs/openms`` for simulation. Run ``sh envs/arda/setup_overlay.sh`` for both
sequencing geometries and ``.venv/bin/python scripts/build_openms_digest_fix.py``
for native MS. Recipe READMEs and hash manifests own exact installation steps.
Do not substitute external repertoire assembly or read mapping for ARDA.

Simulation disk accounting
--------------------------

Simulation monitors count this checkout and the declared workload directory,
including an external Nextflow ``-work-dir``. Unrelated sibling repositories are
excluded; nested workload directories are counted once. CPU, memory and disk
limits remain explicit. This does not run or revalidate a scientific simulation.

Documentation build
-------------------

.. code-block:: sh

   uv sync --extra test --extra docs --frozen
   PATH="$PWD/.venv/bin:$PATH" make -C docs html

The build treats warnings as errors. Open ``docs/_build/html/index.html``.
No scientific analysis runs during the documentation build.
