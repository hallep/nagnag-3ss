''' Convert hg38 chromosome sequence .fasta file into BioSeq record dictionary '''

from lib import ROOT, pickle, SeqIO

# convert to dictionary
hg38_dict = {record.id : record.seq for record in SeqIO.parse(f"{ROOT}/src/hg38.fa", "fasta")}

# save dictionary
pickle.dump(hg38_dict, (open(f"{ROOT}/src/hg38", "wb")))
