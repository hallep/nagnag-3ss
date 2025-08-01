import pandas as pd

CHROMS = list(map(str, range(1,23))) + ["X", "Y", "MT", "M"]

dtypes = {
    "csite_starts" : "str",
    "csite_pos" : "str",
    "csite_inds" : "str",
    "ssite_ind" : "str",
    "nsite_ind" : "str",
    "vsite_ind" : "str",
    "scen_inds" : "str",
    "event_inds" : "str",
    "index" : "str",

    "chrom" : "str",
    "pos" : "int",
    "id" : "str",
    "ref" : "str",
    "alt" : "str",
    "qual" : "str",
    "filter" : "str",
    "info" : "str",

    "off_pos" : "str",
    "snp_pos_create_dbSNP" : "str",
    "snp_pos_alter_dbSNP" : "str",
    "snp_pos_destroy_dbSNP" : "str",
    "snp_pos_create_HGMD_splice" : "str",
    "snp_pos_alter_HGMD_splice" : "str",
    "snp_pos_destroy_HGMD_splice" : "str",
    "snp_pos_create_ClinVar" : "str",
    "snp_pos_alter_ClinVar" : "str",
    "snp_pos_destroy_ClinVar" : "str",

    "snp_id_create_dbSNP" : "str",
    "snp_id_alter_dbSNP" : "str",
    "snp_id_destroy_dbSNP" : "str",
    "snp_id_create_HGMD_splice" : "str",
    "snp_id_alter_HGMD_splice" : "str",
    "snp_id_destroy_HGMD_splice" : "str",
    "snp_id_create_ClinVar" : "str",
    "snp_id_alter_ClinVar" : "str",
    "snp_id_destroy_ClinVar" : "str",

    "snp_id" : "str",
    "vnt_id" : "str",
    "ps_aa_ref" : "str",
    "ps_aa_alt" : "str",
    "ds_aa_ref" : "str",
    "ds_aa_alt" : "str",
}

''' ===== INIT ===== '''

vcf_columns = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info"]

info_fields = {
    "dbSNP" : {
        "RS" : int,
        "RSPOS" : int,
        "RV" : bool,
        "VP" : str,
        "GENEINFO" : str,
        "dbSNPBuildID" : int,
        "SAO" : int,
        "SSR" : int,
        "WGT" : int,
        "VC" : str,
        "PM" : bool,
        "TPA" : bool,
        "PMC" : bool,
        "S3D" : bool,
        "SLO" : bool,
        "NSF" : bool,
        "NSM" : bool,
        "NSN" : bool,
        "REF" : bool,
        "SYN" : bool,
        "U3" : bool,
        "U5" : bool,
        "ASS" : bool,
        "DSS" : bool,
        "INT" : bool,
        "R3" : bool,
        "R5" : bool,
        "OTH" : bool,
        "CFL" : bool,
        "ASP" : bool,
        "MUT" : bool,
        "VLD" : bool,
        "G5A" : bool,
        "G5" : bool,
        "HD" : bool,
        "GNO" : bool,
        "KGPhase1" : bool,
        "KGPhase3" : bool,
        "CDA" : bool,
        "LSD" : bool,
        "MTP" : bool,
        "OM" : bool,
        "NOC" : bool,
        "WTD" : bool,
        "NOV" : bool,
        "CAF" : str,
        "COMMON" : int,
        "TOPMED" : str
    },

    "HGMD" : {
        "CLASS" : str,
        "MUT" : str,
        "GENE" : str,
        "STRAND" : str,
        "DNA" : str,
        "PROT" : str,
        "DB" : str,
        "PHEN" : str,
        "RANKSCORE" : float,
        "SVTYPE" : str,
        "END" : int,
        "SVLEN" : int  
    },
    
    "ClinVar" : {
        "AF_ESP" : float,
        "AF_EXAC" : float,
        "AF_TGP" : float,
        "ALLELEID" : int,
        "CLNDN" : str,
        "CLNDNINCL" : str,
        "CLNDISDB" : str,
        "CLNDISDBINCL" : str,
        "CLNHGVS" : str,
        "CLNREVSTAT" : str,
        "CLNSIG" : str,
        "CLNSIGCONF" : str,
        "CLNSIGINCL" : str,
        "CLNSIGSCV" : str,
        "CLNVC" : str,
        "CLNVCSO" : str,
        "CLNVI" : str,
        "DBVARID" : str,
        "GENEINFO" : str,
        "MC" : str,
        "ONCDN" : str,
        "ONCDNINCL" : str,
        "ONCDISDB" : str,
        "ONCDISDBINCL" : str,
        "ONC" : str,
        "ONCINCL" : str,
        "ONCREVSTAT" : str,
        "ONCSCV" : str,
        "ONCCONF" : str,
        "ORIGIN" : str,
        "RS" : str,
        "SCIDN" : str,
        "SCIDNINCL" : str,
        "SCIDISDB" : str,
        "SCIDISDBINCL" : str,
        "SCIREVSTAT" : str,
        "SCI" : str,
        "SCIINCL" : str,
        "SCISCV" : str
    }
}

info_defaults = {
    int : -1,
    float : -1,
    bool : False,
    str : "."
}

dtypes.update({k:v.__name__ for info in info_fields.values() for k,v in info.items()})

''' ===== FIND ===== '''

''' ----- FILTER ----- '''

# FILTER NAGNAG-CREATING SNPs
def filter_create(df:pd.DataFrame, match_strand:bool) -> pd.DataFrame:

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
    
    if match_strand:
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

# FILTER NAGNAG-ALTERING SNPs
def filter_alter(df:pd.DataFrame, _:bool) -> pd.DataFrame:

    ''' Filter NAGNAG-altering SNPs
    
    Conditions:
    * affected base ("p") is 0 (n1) or 3 (n2)
    '''

    return df[(df["p"] == 0) | (df["p"] == 3)]

# FILTER NAGNAG-DESTROYING SNPs
def filter_destroy(df:pd.DataFrame, _:bool) -> pd.DataFrame:

    ''' Filter NAGNAG-destroying SNPs
    
    Conditions:
    * affected base ("p") is NOT 0 (n1) or 3 (n2)
    '''

    return df[(df["p"] != 0) & (df["p"] != 3)]

''' ----- INFO ----- '''

# NAGNAG Effects ("create", "alter", "destroy"): [0] splice site type, [1] filter, [2] description
effect_info = {
    "create" : {
        0 : "1off",
        1 : filter_create,
        2 : "creating"
    },
    "alter" : {
        0 : "canon",
        1 : filter_alter,
        2 : "altering"
    },
    "destroy" : {
        0 : "canon",
        1 : filter_destroy,
        2 : "destroying"
    }
}
''' NAGNAG variant effects

keys: "create", "alter", "destroy"

values:
* [0] splice site type {str}
* [1] filter {func}
* [2] gerund description {str}
'''

site_info = {
    "1off" : {
        0 : {
            True : "/mnt/data_3/hallep/nagnag/snps/intersect/1off_chr#.bed",
            False : "/mnt/data_3/hallep/nagnag/snps/intersect/1off.bed"
        },
        1 : "/mnt/data_3/hallep/nagnag/snps/1off_splice_sites.txt"
    },
    "canon" : {
        0 : {
            True : "/mnt/data_3/hallep/nagnag/snps/intersect/canon_chr#.bed",
            False : "/mnt/data_3/hallep/nagnag/snps/intersect/canon.bed",
        },
        1 : "/mnt/data_3/hallep/nagnag/nagnag_3ss.txt"
    },
    "all" : {
        0 : {
            True : "/mnt/data_3/hallep/nagnag/snps/intersect/3ss_chr#.bed",
            False : "/mnt/data_3/hallep/nagnag/snps/intersect/3ss.bed"
        },
        1 : "/mnt/data_3/hallep/nagnag/3ss_uppercase.txt"
    }
}
''' Splice site information

keys: "1off", "canon", "all"

values:
    * [0] .bed file {dict[bool, str]}
        * True (by_chrom)
        * False (by_chrom)
    * [1] .txt file
'''

var_info = {
    "dbSNP" : {
        0 : True,
        1 : "dbsnp_chr#.bed",
        2 : "/mnt/data_3/hallep/reference/snp_data/dbSNP/chrom_split/dbSNP_chr#.txt",
        3 : "id",
        4 : False,
        5: info_fields["dbSNP"]
    },
    "dbSNP_common" : {
        0 : True,
        1 : "dbsnp_chr#.bed",
        2 : "/mnt/data_3/hallep/reference/snp_data/dbSNP/chrom_split/dbSNP_chr#.txt",
        3 : "id",
        4 : False,
        5: info_fields["dbSNP"]
    },
    "dbSNP_rare" : {
        0 : True,
        1 : "dbsnp_chr#.bed",
        2 : "/mnt/data_3/hallep/reference/snp_data/dbSNP/chrom_split/dbSNP_chr#.txt",
        3 : "id",
        4 : False,
        5: info_fields["dbSNP"]
    },
    "ClinVar" : {
        0 : False,
        1 : "clinvar.bed",
        2 : "/mnt/data_3/hallep/reference/snp_data/ClinVar/canon_clinvar_snps.txt",
        3 : "id",
        4 : False,
        5 : info_fields["ClinVar"]
    },
    "HGMD" : {
        0 : False,
        1 : "hgmd_vcf.bed",
        2 : "/mnt/data_3/hallep/reference/snp_data/HGMD/all_HGMD_snps.txt",
        3 : "id",
        4 : False,
        5 : info_fields["HGMD"]
    },
    "HGMD_splice" : {
        0 : False,
        1 : "hgmd_splice.bed",
        2 : "/mnt/data_3/hallep/reference/snp_data/HGMD/hgmd_splice_variants.txt",
        3 : "acc_num",
        4 : True,
        5 : None
    },
}
''' Variant information

keys: "dbSNP", "dbSNP_common", "dbSNP_rare", ClinVar", "HGMD", "HGMD_splice"

values:
* [0] by_chrom {bool}
* [1] bed {str}
* [2] src_txt {str}
* [3] info_col {str}
* [4] match_strand {bool}
* [5] info_fields {dict}
'''

''' ===== ANALYZE ===== '''

snp_file = "/mnt/data_3/hallep/nagnag/snps/[SNP]_[SITE]_snps.txt"
site_file = "/mnt/data_3/hallep/nagnag/snps/[SNP]_[SITE]_snp_containing_ssites.txt"

snp_dbs = ["dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"]
site_types = ["3ss", "nagnag"]

''' ===== SCENARIOS ===== '''

coord_col = {
    "dbSNP_common" : "pos",
    "dbSNP_rare" : "pos",
    "ClinVar" : "pos",
    "HGMD_splice" : "start"
}

maf_col = {

    "dbSNP_common" : {
        "src" : ["COMMON", "G5", "G5A", "maf_TGP", "maf_TOPMed"],
        "dst" : ["COMMON", "G5", "G5A", "maf_TGP", "maf_TOPMed"]
    },

    "dbSNP_rare" : {
        "src" : ["COMMON", "G5", "G5A", "maf_TGP", "maf_TOPMed"],
        "dst" : ["COMMON", "G5", "G5A", "maf_TGP", "maf_TOPMed"]
    },

    "ClinVar" : {
        "src" : ["AF_ESP", "AF_EXAC", "AF_TGP"],
        "dst" : ["maf_ESP", "maf_ExAC", "maf_TGP"]
    },

    "HGMD_splice" : {
        "src" : [],
        "dst" : []
    }
}

nn_effects = ["create", "alter_ref", "alter_alt", "destroy"]
