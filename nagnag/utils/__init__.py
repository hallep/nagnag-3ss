''' Utilities
* `ROOT`: root filepath
* `seq`: genomic sequence variables and functions
* `vnt`: genomic variant variables and functions
'''

__all__ = ["ROOT", "seq", "vnt", "log_script", "log_fn"]

# Root Folder
from pathlib import Path
ROOT = Path(__file__).parent.parent.parent

# Modules
from . import seq
from . import vnt

# Run Log
from .log import log_script, log_fn
