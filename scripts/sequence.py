''' Helper functions pertaining to DNA, RNA, and amino acid sequences '''

from lib import ROOT, pickle, Seq

hg38 = pickle.load(open(f"{ROOT}/src/hg38", "rb"))

def ss_seq(chrom:str, pos:int, strand:str, up:int, down:int, upper:bool=False) -> str:

    ''' Get splice site sequence:
    
    {up} upstream nucleotides + {down} downstream nucleotides
    '''

    seq = ""

    # Plus (+) Strand
    if strand == "+":
        seq = str(hg38[chrom][pos-up:pos+down])

    # Minus (-) Strand
    elif strand == "-":
        seq = str(hg38[chrom][pos-down:pos+up].reverse_complement())

    # convert to uppercase
    if upper:
        return seq.upper()

    return seq
