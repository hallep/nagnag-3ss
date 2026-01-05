''' Utilities
* `ROOT`: root directory filepath
* `seq`: genomic sequence variables and functions
* `vnt`: genomic variant variables and functions
* `log_script`, `log_fn`: debug log functions
* `single_map`, `multi_map`: parallel computing functions
'''

__all__ = ["ROOT", "seq", "vnt", "log_script", "log_fn", "single_map", "multi_map"]

# Root Folder
from pathlib import Path
ROOT = Path(__file__).parent.parent.parent
''' Root directory '''

# Modules
from . import seq
from . import vnt

# Run Log
from .log import log_script, log_fn

# parallel processing
from .compute import single_map, multi_map
