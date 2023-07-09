# crnpy/__init__.py
""" `crnpy` is a package for processing observations from cosmic ray neutron detectors.

Modules exported by this package:

- `crnpy`: Provide several functions to process observations from cosmic ray neutron detectors.
"""

# GEt version from setup.py
from pathlib import Path
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "unknown"

# Import all functions from crnpy

from .crnpy import *

