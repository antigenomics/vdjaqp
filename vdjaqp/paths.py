"""Portable data roots shared with the ``vdjaqp-data`` bootstrap.

``VDJAQP_DATA`` defaults to ``./data``. The external cohort layout is under
``$VDJAQP_DATA/external/gut``; see the dataset reference for access and provenance.
"""
import os
from pathlib import Path


def data_root():
    return Path(os.environ.get('VDJAQP_DATA', 'data')).expanduser()


DATA = data_root()
EXTERNAL = DATA / 'external'
GUT = EXTERNAL / 'gut'
