"""Report the installed scoring runtime without loading private inputs or models."""
from importlib import import_module
from importlib.metadata import version
import sys

for module in ('vdjaqp.pipeline', 'vdjaqp.match', 'vdjaqp.pgen', 'vdjaqp.corpus'):
    import_module(module)
print('Python', sys.version.split()[0])
for package in ('vdjaqp', 'polars', 'numpy', 'pyarrow', 'seqtree', 'arda-mapper', 'vdjtools'):
    print(package, version(package))
