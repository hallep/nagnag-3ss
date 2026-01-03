''' Helper functions for variant processing '''

from . import ROOT
from lib import pd

dtypes = {

    # 1off_splice_sites.txt
    "csite_starts" : "str",
    "csite_pos" : "str",
    "csite_inds" : "str",

    # vcf columns
    "chrom" : "str",
    "pos" : "int",
    "id" : "str",
    "ref" : "str",
    "alt" : "str",
    "qual" : "str",
    "filter" : "str",
    "info" : "str",


}

vcf_cols = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info"]
txt_cols = ["id", "chrom", "pos", "ref", "alt", "qual", "filter", "info"]

# --- Variant Databases --- #
class VntDatabase:
    def __init__(self, db:str, by_chrom:str, stranded:bool, split_freq:bool):
        self.by_chrom = by_chrom
        self.stranded = stranded
        self.split_freq = split_freq
        self.txt = f"{ROOT}/variants/src/{db}{"_$" if split_freq else ""}{"_#" if by_chrom else ""}.txt"
        self.bed = f"{ROOT}/variants/bed/{db}{"_$" if split_freq else ""}{"_#" if by_chrom else ""}.bed"

DB = {
    "dbSNP" : VntDatabase("dbSNP", True, False, True),
    "dbSNP_common" : VntDatabase("dbSNP_common", True, False, False),
    "dbSNP_rare" : VntDatabase("dbSNP_rare", True, False, False),
    "ClinVar" : VntDatabase("ClinVar", False, False, False),
    "HGMD_splice" : VntDatabase("HGMD_splice", False, True, False),
}
''' Database Information
* keys: "dbSNP", "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"
* values: VntDatabase
    * by_chrom (bool): whether to search by chromosome
    * stranded (bool): whether strand matters
    * split_freq (bool): whether to split by COMMON flag
    * txt (str): path to .txt file
    * bed (str): path to .bed file
'''

all_vnt_dbs = ["dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"]
''' List of all variant databases: dbSNP_common, dbSNP_rare, ClinVar, HGMD_splice '''

# --- Splice Sites --- #
class SpliceVntSite:
    def __init__(self, txt:str, bed:str):
        self.txt = f"{ROOT}/sites/{txt}.txt"
        self.bed = {
            True : f"{ROOT}/variants/bed/{bed}_#.bed",
            False : f"{ROOT}/variants/bed/{bed}.bed"
        }

SS = {
    "1off" : SpliceVntSite("1off_splice_sites", "1off"),
    "canon" : SpliceVntSite("nagnag_3ss", "canon"),
    "3ss" : SpliceVntSite("3ss_uppercase", "3ss")
}
''' Splice site information
* keys: "1off", "canon", "all"
* values: SpliceVntSites
    * txt (str): path to .txt file
    * bed (str): dictionary of paths to .bed files (keys are True/False with respect to by_chrom)
'''

# --- NAGNAG Effects --- #
def filter_effect(df:pd.DataFrame, effect:str=None, stranded:bool=None) -> pd.DataFrame:

    ''' Filter found'''
    # no filtering
    if effect == None:
        return df

    # create
    def filter_create(df:pd.DataFrame, stranded:bool) -> pd.DataFrame:

        ''' Filter NAGNAG-creating SNPs
        
        Conditions:
        * match_strand = True:
            * off-position ("p") is 1 (a1) or 4 (a2)
                    * strand ("r") is "+": alt allele ("a") must be "A"
                    * strand ("r") is "-": alt allele ("a") must be "T"
            * off-position ("p") is 2 (g1) or 5 (g2)
                    * strand ("r") is "+": alt allele ("a") must be "G"
                    * strand ("r") is "-": alt allele ("a") must be "C"
        '''
        
        if stranded:
            b = ((df["p"] == 1) | (df["p"] == 4)) & (df["a"].str.split(",").apply(lambda x : "A" in x))
            h = ((df["p"] == 2) | (df["p"] == 5)) & (df["a"].str.split(",").apply(lambda x : "G" in x))

            return df[b | h]

        # off position is "A"
        bp = ((df["p"] == 1) | (df["p"] == 4)) & (df["r"] == "+") & (df["a"].str.split(",").apply(lambda x : "A" in x))
        bm = ((df["p"] == 1) | (df["p"] == 4)) & (df["r"] == "-") & (df["a"].str.split(",").apply(lambda x : "T" in x))

        # off position is "G"
        hp = ((df["p"] == 2) | (df["p"] == 5)) & (df["r"] == "+") & (df["a"].str.split(",").apply(lambda x : "G" in x))
        hm = ((df["p"] == 2) | (df["p"] == 5)) & (df["r"] == "-") & (df["a"].str.split(",").apply(lambda x : "C" in x))
        
        return df[bp | bm | hp | hm]

    # alter
    def filter_alter(df:pd.DataFrame, _:bool) -> pd.DataFrame:

        ''' Filter NAGNAG-altering SNPs
        
        Conditions:
        * affected base ("p") is 0 (n1) or 3 (n2)
        * has at least one valid (i.e., non ".") base
        '''

        v = df["a"].str.split(",").apply(lambda x: len(set(x).intersection(["A","C","G","T"])) > 0)
        return df[((df["p"] == 0) | (df["p"] == 3)) & (v)]

    # destroy
    def filter_destroy(df:pd.DataFrame, _:bool) -> pd.DataFrame:

        ''' Filter NAGNAG-destroying SNPs
        
        Conditions:
        * affected base ("p") is NOT 0 (n1) or 3 (n2)
        '''

        return df[(df["p"] != 0) & (df["p"] != 3)]

    # filter function
    fn = {
        "create" : filter_create,
        "alter" : filter_alter,
        "destroy" : filter_destroy,
    }[effect]

    return fn(df, stranded)

class VntEffect:
    def __init__(self, site_type:str, adj:str):
        self.site_type = site_type
        self.adj = adj

EFF = {
    "create" : VntEffect("1off", "creating"),
    "alter" : VntEffect("canon", "altering"),
    "destroy" : VntEffect("canon", "destroying"),
}
''' NAGNAG variant effects 
* keys: "create", "alter", "destroy"
* values: VntDatabase
    * site_type (str): "1off" or "canon"
    * adj (str): adjective form of effect
'''
