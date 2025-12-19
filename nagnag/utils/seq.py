''' Helper functions pertaining to DNA, RNA, and amino acid sequences '''

from . import ROOT
from lib import pickle, re, Seq

hg38 = pickle.load(open(f"{ROOT}/src/hg38", "rb"))

# ===== Sequence ===== #

CHROMS = [f"chr{c}" for c in list(map(str, range(1, 23))) + ["X", "Y", "M"]]
N = ["A", "C", "G", "T"]
AA = ["A", "F", "I", "L", "M", "P", "V", "W", "C", "N", "Q", "S", "T", "Y", "D", "E", "H", "K", "R", "G", "*"]

# Splice Site Sequence
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

# Cleavage Site Sequence
def get_cs_seq(chrom:str, pos:int, strand:str, upper:bool=False) -> str:

    ''' Get cleavage site sequence (prececeding 3 nucleotides) '''

    return get_ss_seq(chrom, pos, strand, 3, 0, upper)

# ===== Sites ===== #

def is_nag(seq:str) -> bool:
    return re.fullmatch(".AG", seq, re.IGNORECASE) != None

# Splice Type
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

# ===== Proteome Effects ===== #

# Transition
TSNS = ["vc_ps", "vc_ds", "aa_ps", "aa_ds", "aattype"]
AATTYPES = {
    "DID" : "Duplicate Indel",
    "NID" : "Pre-Indel",
    "CID" : "Post-Inde",
    "IDR" : "Indel + Replacement",
    "NC" : "No Change",
    "ET" : "Elongation/Truncation"
}

def get_aattype(ps_aa:str, ds_aa:str) -> str:

    ''' Amino acid transition type '''

    # non-coding
    if (ps_aa == "_") and (ds_aa == "_"):
        return "_"
    
    # phase 0
    if ds_aa == ".":
        return ps_aa

    # phases 1 and 2
    pn, pc = ps_aa
    d = ds_aa
        
    # NC (no change)
    if (pn == "*") and (d == "*"):
        return "NC"

    # ET (elongation/truncation)
    if (((pn != "*") and (pc != "*") and (d == "*")) or 
        (((pn == "*") or (pc == "*")) and (d != "*"))):
        return "ET"
    
    # DID (duplicate indel)
    if (pn == d) and (pc == d):
        return "DID"

    # NID (pre indel)
    if (pn != d) and (pc == d):
        return "NID"
    
    # CID (post indel)
    if (pn == d) and (pc != d):
        return "CID"
    
    # IDR (indel + replacement)
    return "IDR"

def get_tsn(up:str, motif:str, down:str, phase:str|int) -> tuple[str, str, str, str, str]:

    ''' Get variable RNA codons and amino acids '''

    # non-coding
    if str(phase) == "-1":
        ps_codon = "_"
        ps_aa = "_"

        ds_codon = "_"
        ds_aa = "_"

    # phase 0
    elif str(phase) == "0":
        ps_codon = motif[3:6]
        ps_aa = str(Seq.translate(ps_codon))

        ds_codon = "."
        ds_aa = "."

    # phase 1
    elif str(phase) == "1":
        ps_codon = up[-1] + motif[3:6] + down[:2]
        ps_aa = str(Seq.translate(ps_codon))

        ds_codon = up[-1] + down[:2]
        ds_aa = str(Seq.translate(ds_codon))

    # phase 2
    elif str(phase) == "2":
        ps_codon = up[-2:] + motif[3:6] + down[0]
        ps_aa = str(Seq.translate(ps_codon))

        ds_codon = up[-2:] + down[0]
        ds_aa = str(Seq.translate(ds_codon))
    
    aattype = get_aattype(ps_aa, ds_aa)

    return ps_codon, ds_codon, ps_aa, ds_aa, aattype

# Amino Acid Outcomes
AA = ["A", "F", "I", "L", "M", "P", "V", "W", "C", "N", "Q", "S", "T", "Y", "D", "E", "H", "K", "R", "G", "*"]
OUT = ["ins", "del"]

def inserted(p:int, ps:str, ds:str) -> str:

    ''' Insertion: which amino acid(s) (incl. *) are inserted (from distal to proximal)
    * if duplicate indel: return {ds}
    * if pre-indel, return A_PN
    * if post-indel, return A_PC
    * if indel + replcement, return {ps}
    * if elongation, return {ps}
    * if truncation, return "*"
    * if no insertion (no change), return "."
    '''

    # non-coding
    if p == -1:
        return "_"
    
    # phase 0
    if p == 0:
        return ps
    
    if ps[0] == "*":
        ps = ps[0]
    
    # no change
    if ps == ds:
        return "."

    # indel
    aa_ps = list(ps)
    
    if ds in aa_ps:
        aa_ps.remove(ds)

    return "".join(aa_ps)

def deleted(p:int, ps:str, ds:str) -> str:

    ''' Deletion: which amino acid (incl. *) is deleted (from distal to proximal)
    * if deletion, return {ds}
    * if no deletion, return "."
    '''

    # non-coding
    if p == -1:
        return "_"

    # phase 0
    if p == 0:
        return "."

    # indel, no replacement
    if ds in ps:
        return "."
    
    # indel + replacement
    return ds

def get_categorical_outcomes(phase:int, ps_aa:str, ds_aa:str) -> tuple[str, str]:

    ''' Get categorical amino acid outcomes (inserted, deleted) '''

    # non-coding
    if phase == -1:
        return "_", "_"

    # categorical outcomes
    insert = inserted(phase, ps_aa, ds_aa)
    delete = deleted(phase, ps_aa, ds_aa)

    return insert, delete
