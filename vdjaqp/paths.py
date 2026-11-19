"""Local data roots for the manuscript analyses (``analyses/``).

``VDJAQP_DATA`` (default ``~/work/academy/2026/vdjaqp``) is the local workspace; the
public gut cohorts used in the paper sit in ``$VDJAQP_DATA/external/gut`` with the
layout of the former ``colitis_igprotseq`` workspace (``delivery/``, ``raw/``,
``references/``, ``outcome_inventory/``). Processed files are also in the HF dataset
``isalgo/igprotseq`` (``vdjaqp-data --group external``). Historical receipts keep the
original absolute paths in the private manuscript records.
"""
import os
from pathlib import Path

DATA = Path(os.environ.get('VDJAQP_DATA', '~/work/academy/2026/vdjaqp')).expanduser()
EXTERNAL = DATA / 'external'
GUT = EXTERNAL / 'gut'
