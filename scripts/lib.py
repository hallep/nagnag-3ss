__all__ = ["ROOT", "requests", "io", "gzip", "pickle", "SeqIO"]

from pathlib import Path
ROOT = Path(__file__).parent.parent

# Input/Output
import requests
import io
import gzip
import pickle
from Bio import SeqIO

