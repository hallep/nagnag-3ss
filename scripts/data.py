''' Convert hg38 chromosome sequence .fasta file into BioSeq record dictionary;

**Created Files:**
* `src/hg38`
'''

from utils import ROOT
from utils.lib import pickle, SeqIO

print("parsing GRCh38/hg38 .fasta file...", end="", flush=True)

# convert to dictionary
hg38_dict = {record.id : record.seq for record in SeqIO.parse(f"{ROOT}/src/hg38.fa", "fasta")}

# save dictionary
pickle.dump(hg38_dict, (open(f"{ROOT}/src/hg38", "wb")))

print("done")
