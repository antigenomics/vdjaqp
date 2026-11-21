# Workflow entry points

The entry point is `../main.nf`. See the [pipeline guide](../docs/pipeline.rst),
[parameters](../docs/parameters.rst), [input/output contracts](../docs/input-output.rst)
and [input preparation](../docs/preparation.rst).

Copy `samples.example.csv` and supply existing input tables and a fresh output
directory. The `analysis_node` and `analysis_slurm` profiles bound concurrent
analysis by explicit CPU and memory budgets. Site accounts and storage paths
belong in local configuration.

Scientific acceptance records, manuscript scheduler receipts and recovery scripts
are retained privately with the manuscript. The software tests use invented data.

Snakemake also accepts the observed-analysis sample sheet. See the [pipeline
guide](../docs/pipeline.rst) for the command and Slurm settings.
