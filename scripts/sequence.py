''' Helper functions pertaining to DNA, RNA, and amino acid sequences '''

from lib import ROOT, pickle, re, Seq

hg38 = pickle.load(open(f"{ROOT}/src/hg38", "rb"))

def get_ss_seq(chrom:str, pos:int, strand:str, up:int=0, down:int=0, upper:bool=False) -> str:

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

def get_cs_seq(chrom:str, pos:int, strand:str, upper:bool=False) -> str:

    ''' Get cleavage site sequence (prececeding 3 nucleotides) '''

    return get_ss_seq(chrom, pos, strand, 3, 0, upper)

def is_nag(seq:str) -> bool:
    return re.fullmatch(".AG", seq, re.IGNORECASE) != None

def get_sstype(seq:str) -> str:

    ''' Get splice site type (based on number and canonicity of NAGs) '''

    # canonical 1-NAG (1C)
    if re.fullmatch(".AG", seq, re.IGNORECASE):
        return "1C"
    
    # non-canonical 1-NAG (1NC)
    if len(seq) == 3:
        return "1NC"

    # canonical NAGNAG (2C)
    if re.fullmatch(".AG.AG", seq, re.IGNORECASE):
        return "2C"
    
    # non-canonical 2-NAG (2C)
    if len(seq) == 6:
        return "2NC"
    
    # 3+ NAG
    if len(seq) % 3 == 0:
        return "3+"

    return "."
