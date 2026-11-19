"""Link LC-MS/MS antibody peptides to BCR repertoire clonotypes."""

from .resources import DEFAULT_CPU_BUDGET, cpu_budget, resolve_cpu_budget

resolve_cpu_budget()

__version__ = "0.1.0"
__all__ = ["DEFAULT_CPU_BUDGET", "cpu_budget", "resolve_cpu_budget", "__version__"]
