''' Helper functions pertaining to DNA and amino acid sequences '''

import pickle
from Bio import Seq

# Get sequence dictionary
with open("source/hg38", "rb") as file:
    hg38 = pickle.load(file)

# Get splice site sequences
def get_ss_seq(chrom:str, pos:int, strand:str, up:int, down:int, upper=False) -> str:

    ''' Get sequences surrounding a cleavage site
    
    {up} upstream nucleotides + {down} downstream nucleotides
    
    Parameters
    ----------
    chrom : str
        chromosome name: chr#
    pos : int
        0-indexed chromosome position of cleavage site
         * plus (+) strand:
             * 3': pos is first base of downstream exon
             * 5': pos is first base of downstream intron
         * minus (-) strand:
             * 3': pos is last base of upstream intron
             * 5': pos is last base of upstream exon
    strand : str {"+", "-"}
        strand on which cleavage site is found
         * "+": get direct sequence
         * "-": get reverse complement of sequence
    up : int
        number of upstream nucleotides to get
         * 3': nucleotides are intronic
         * 5': nucleotides are exonic
    down : int
        number of downstream nucleotides to get
         * 3': nucleotides are exonic
         * 5': nucleotides are intronic
    upper : bool (default = False)
        whether to ignore repeat masker (lowercase)
         * True: convert sequence to all uppercase
         * False: leave sequence as is

    Returns
    -------
    tri : str
        the nucleotides (on the specified strand) surrounding a cleavage site
    '''
    
    tri = ""

    # plus (+) strand
    if strand == "+":
        tri = str(hg38[chrom][pos-up:pos+down])

    # minus (-) strand    
    elif strand == "-":
        tri = str(hg38[chrom][pos-down:pos+up].reverse_complement())
    
    # convert to uppercase
    if upper:
        return tri.upper()
    
    return tri

# Get trinucleotide preceeding 3' splice site 
def get_3cs_tri(chrom:str, pos:int, strand:str, upper=False) -> str:

    ''' Get the intronic trinucleotide sequence preceeding a 3' cleavage site
    
    Call get_ss_seq(chrom, pos, strand, 3, 0, upper)

    Parameters
    ----------
    chrom : str
        chromosome name: chr#
    pos : int
        0-indexed chromosome position of cleavage site
         * plus (+) strand: pos is first base of downstream exon
         * minus (-) strand: pos is last base of upstream intron
    strand : str {"+", "-"}
        strand on which cleavage site is found
         * "+": get direct sequence
         * "-": get reverse complement of sequence
    upper : bool (default = False)
        whether to ignore repeat masker (lowercase)
         * True: convert sequence to all uppercase
         * False: leave sequence as is

    Returns
    -------
    _ : str
        the 3 nucleotides (on the proper strand) preceeding a cleavage site
    '''

    return get_ss_seq(chrom=chrom, pos=pos, strand=strand, up=3, down=0, upper=upper)

# Get NAGNAG motifs
def get_NAGNAGs() -> list:
    
    ''' Get list of NAGNAG motifs
    
    Returns
    -------
    _ : list
        all NAGNAG motifs, with N1 = [A, C, G, T], then N2 = [A, C, G, T]
    '''

    n = ["A", "C", "G", "T"]

    return [f"{n1}AG{n2}AG" for n1 in n for n2 in n]

# Calculate amino acid transition
def get_aa_transition(up:str, n2:str, down:str, frame:str|int) -> tuple[str, str, str, str]:

    ''' Get amino acid transition of a NAGNAG
    
    Parameters
    ----------
    up : str
        last 3 bases of upstream exon
    n2 : str {"A", "C", "G", "T"}
        identity of N2
    down : str
        first 3 bases of downstream exon
    frame : str {"-1", "0", "1", "2"} or int {-1, 0, 1, 2}
        frame of splice scenario

    Returns
    -------
    ps_codon : str
        variable codons when proximally spliced
         * phase -1: "_"
         * phase 0: 3 bases (N2 A G)
         * phase 1: 6 bases (U-1 N2 A G D1 D2)
         * phase 2: 6 bases (U-2 U-1 N A G D1)
    ps_aa : str
        variable amino acid(s) when proximally spliced
         * phase -1: "_"
         * phase 0: 1 amino acid
         * phase 1/2: 2 amino acids
    ds_codon : str
        variable codons when distally spliced
         * phase -1: "_"
         * phase 0: "."
         * phase 1: 3 bases (U-1 D1 D2)
         * phase 2: 3 bases (U-2 U-1 D1)
    ds_aa : str
        variable amino acid when distally spliced
         * phase -1: "_"
         * phase 0: "."
         * phase 1/2: 1 amino acid
    '''

    # phase -1 (not CDS)
    if str(frame) == "-1":
        ps_codon = "_"
        ps_aa = "_"

        ds_codon = "_"
        ds_aa = "_"

    # phase 0
    elif str(frame) == "0":
        ps_codon = n2 + "AG"
        ps_aa = str(Seq.translate(ps_codon))

        ds_codon = "."
        ds_aa = "."

    # phase 1
    elif str(frame) == "1":
        ps_codon = up[-1] + n2 + "AG" + down[:2]
        ps_aa = str(Seq.translate(ps_codon))

        ds_codon = up[-1] + down[:2]
        ds_aa = str(Seq.translate(ds_codon))

    # phase 2
    elif str(frame) == "2":
        ps_codon = up[-2:] + n2 + "AG" + down[0]
        ps_aa = str(Seq.translate(ps_codon))

        ds_codon = up[-2:] + down[0]
        ds_aa = str(Seq.translate(ds_codon))
    
    return ps_codon, ps_aa, ds_codon, ds_aa
