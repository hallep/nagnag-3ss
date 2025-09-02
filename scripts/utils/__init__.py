''' Run utilities
* `ROOT`: root filepath
* `lib`: external Python packages
* `seq`: genomic sequence variables and functions
* `vnt`: genomic variant variables and functions
'''

__all__ = ["ROOT", "lib", "seq", "vnt", "run"]

# Root Folder
from pathlib import Path
ROOT = Path(__file__).parent.parent.parent

from . import seq
from . import vnt

def run(fn, *args, desc:str, sub:bool=False):
    print(f"{" - " if sub else ""}{desc}...", end="", flush=True)
    fn(*args)
    print("done")
