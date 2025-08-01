import pandas as pd
import numpy as np
from parallelbar import progress_map, progress_imap, progress_starmap
from snps_init import get_info_value
from snps_tools import CHROMS, dtypes, effect_info, site_info, var_info, snp_file, site_file, snp_dbs, site_types

''' ===== SNP Effects ===== '''

# NAGNAG EFFECT
def get_nagnag_effect(snp_txt:pd.DataFrame, site_txt:pd.DataFrame, stranded:bool):

    ''' Get SNP NAGNAG effect (create, alter, destroy)
    
    added columns:
    * "snp_type" : str {"CREATE", "ALTER", "DESTROY"}
        SNP effect in terms of NAGNAGs
    * "alt_ind" : int
        0-index of allele (in "alt") for SNP type
        * if ALTER or DESTROY, use 0 (first alt allele)
        * if CREATE, must be "A" (pos = 1 or 4) or "G" (pos = 2 or 5)
    * "ref_single", "alt_single" : int {0, 1}
        whether the reference and alternate alleles (respectively) are a single base or not

    Parameters
    ----------
    snp_txt : str
        filepath to .txt file of SNPs
    site_txt : str
        filepath to .txt file of splice sites
    stranded : bool
        whether strand was considered when finding intersecting SNPs
    '''

    # load SNPs
    snps = pd.read_csv(snp_txt, sep="\t", index_col=0, dtype=dtypes)

    # load splice sites
    sites = pd.read_csv(site_txt, sep="\t", index_col=0, dtype=dtypes)

    types = [set() for _ in snps.index]
    alt_ind = [-1] * len(snps)
    single = [0] * len(snps)

    # for each SNP:
    for x,(snp,alts,inds) in enumerate(zip(snps.index.values, snps["alt"].values, snps["vsite_ind"].values)):

        info = set()
        
        # alternative alleles
        alts = alts.upper().split(",")
        alt_lens = {len(a) : ai for ai,a in enumerate(alts[::-1]) if a != "."}

        # for all affected sites:
        for i in str(inds).split(","):

            site = sites.loc[i]

            p,t = {n : (int(p),t) for n,p,t in zip(site["snp_ids"].split(","), site["snp_pos"].split(","), site["snp_types"].split(","))}[snp]

            # NAGNAG-altering/destroying
            if t != "CREATE":
                
                if 1 in alt_lens:
                    a = alt_lens[1]
                else:
                    a = 0
                    while alts[a] == ".":
                        a += 1

                        if a >= len(alts):
                            a = -1
                            break

                s = 0 if a == -1 else int(len(alts[a]) == 1)
            
            # NAGNAG-creating
            else:
                
                # create[stranded][site strand][snp pos] = required alt allele for NAGNAG creation
                create = {

                    # stranded == True
                    True : {
                        
                        # site on plus strand
                        "+" : {
                            1 : "A",
                            2 : "G",
                            4 : "A",
                            5 : "G"
                        },

                        # site on minus strand
                        "-" : {
                            1 : "A",
                            2 : "G",
                            4 : "A",
                            5 : "G"
                        }
                    },

                    # stranded == False
                    False : {

                        # site on plus strand
                        "+" : {
                            1 : "A",
                            2 : "G",
                            4 : "A",
                            5 : "G"
                        },

                        # site on minus strand
                        "-" : {
                            1 : "T",
                            2 : "C",
                            4 : "T",
                            5 : "C"
                        }
                    }
                }

                a = alts.index(create[stranded][site["strand"]][p])
                s = 1
            
            # alt allele info
            info.add((t, a, s))
        
        # info
        t, a, s = zip(*info)
        types[x] = ",".join(t)
        alt_ind[x] = ",".join(map(str, a))
        single[x] = ",".join(map(str, s))

    # add columns to DataFrame    
    snps["snp_type"] = types
    snps["alt_ind"] = alt_ind
    snps["ref_single"] = (snps["ref"].str.len() == 1).astype(int)
    snps["alt_single"] = single
    snps.to_csv(snp_txt, sep="\t", index=True)

    print(snps)

# GENERALIZED FUNCTION
def nagnag_effect(db:str):

    ''' Determine NAGNAG-affecting SNPs
    
    call get_nagnag_effect()

    Parameters
    ----------
    db : str {"dbSNP", "ClinVar", "HGMD_splice}
        SNP database
    '''

    get_nagnag_effect(snp_txt=f"/mnt/data_3/hallep/nagnag/snps/found/{db}_nagnag_snps.txt", 
                      site_txt=f"/mnt/data_3/hallep/nagnag/snps/found/{db}_nagnag_snp_containing_ssites.txt",
                      stranded=var_info[db][4])

''' ===== SINGLE BASE SUBSTITUTIONS ===== '''

# SINGLE 3' SPLICE SITE SNPs
def single_3ss_snps(snps:pd.DataFrame, sites:pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:

    ''' Get single base substitution 3' splice site SNPs

    Parameters
    ----------
    snps : pandas.DataFrame
        3'ss-affecting SNPs \\
        required columns: "ref", "alt"
    sites : pandas.DataFrame
        3'ss SNP-containing sites \\
        required column: "snp_id_{db}" for all SNP databases
    
    Returns
    -------
    _ : pandas.DataFrame
        3'ss-affecting SNPs that are single-base substitutions
    _ : pandas.DataFrame
        sites that contain single-base substitution 3'ss-affecting SNPs
    '''

    # SNPs
    single_snps = snps[(snps["ref"].str.len() == 1) & (snps["alt"].str.split(",").apply(lambda alts: any(len(a) == 1 for a in alts)))]

    # sites
    conditions = pd.Series(False, index=sites.index)
    indices = set(single_snps.index.values)
    for col in [f"snp_id_{db}" for db in ["dbSNP", "ClinVar", "HGMD_splice"] if f"snp_id_{db}" in sites.columns]:
        conditions |= sites[col].str.split(",").apply(lambda x: any(i in indices for i in x))
    single_sites = sites[conditions]
    
    return single_snps, single_sites

# SINGLE NAGNAG SNPs
def single_nagnag_snps(snps:pd.DataFrame, sites:pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:

    ''' Get single base substitution NAGNAG SNPs

    Parameters
    ----------
    snps : pandas.DataFrame
        NAGNAG-affecting SNPs \\
        required columns: "ref_single", "alt_single", "snp_type"
    sites : pandas.DataFrame
        NAGNAG SNP-containing sites \\
        required column: "snp_ids", "alt_single", "snp_type"
    
    Returns
    -------
    _ : pandas.DataFrame
        NAGNAG-affecting SNPs that are single-base substitutions
    _ : pandas.DataFrame
        sites that contain single-base substitution NAGNAG-affecting SNPs
    '''

    # SNPs
    single_snps = snps[(snps["ref_single"].astype(str) == "1") & (snps["alt_single"].astype(str) == "1") & (snps["snp_type"].str.contains(",") == False)]

    # sites
    single_sites = sites[sites["snp_ids"].apply(lambda ids: sum([int(i in single_snps.index) for i in ids.split(",")]) > 0)]

    return single_snps, single_sites

# SINGLE BASE SUBSTITUTION SNPs
def isolate_single_snps(src_snps:str, src_sites:str, dst_snps:str, dst_sites:str, fn):

    ''' Isolate single base substitution SNPs

    Parameters
    ----------
    src_snps, src_sites : str
        filepath to source .txt files of all SNPs and SNP-containing sites (respectively)
    dst_snps, dst_sites : str
        filepath to destination .txt file of single SNPs and single SNP-containing sites (respectively)
    fn : function
        isolates single SNPs \\
        * parameter: database (str), all snps (pandas.DataFrame), all sites (pandas.DataFrame)
        * returns: single snps (pandas.DataFrame) and sites (pandas.DataFrame)
    '''
    
    # load SNPs and sites
    all_snps = pd.read_csv(src_snps, sep="\t", index_col=0, dtype=dtypes)
    all_sites = pd.read_csv(src_sites, sep="\t", index_col=0, dtype=dtypes)

    # isolate single SNPs
    single_snps, single_sites = fn(all_snps, all_sites)

    # save
    single_snps.to_csv(dst_snps, sep="\t", index=True)
    single_sites.to_csv(dst_sites, sep="\t", index=True)

    print("\n===== SNPs =====")
    print(single_snps)

    print("\n===== SNPs =====")
    print(single_sites)

# ISOLATE SINGLE 3' SPLICE SITE SNPs: GENERALIZED
def isolate_single_3ss_snps(db:str):

    ''' Isolate single base substitution 3' splice site-affecting SNPs
    
    Parameters
    ----------
    db : str {"dbSNP", "ClinVar", "HGMD_splice}
        SNP database
    '''

    isolate_single_snps(src_snps=f"/mnt/data_3/hallep/nagnag/snps/found/{db}_3ss_snps.txt",
                        src_sites=f"/mnt/data_3/hallep/nagnag/snps/found/{db}_3ss_snp_containing_ssites.txt",
                        dst_snps=f"/mnt/data_3/hallep/nagnag/snps/{db}_3ss_snps.txt",
                        dst_sites=f"/mnt/data_3/hallep/nagnag/snps/{db}_3ss_snp_containing_ssites.txt",
                        fn=single_3ss_snps)

# ISOLATE SINGLE NAGNAG SNPs: GENERALIZED
def isolate_single_nagnag_snps(db:str):

    ''' Isolate single base substitution NAGNAG-affecting SNPs
    
    Parameters
    ----------
    db : str {"dbSNP", "ClinVar", "HGMD_splice}
        SNP database
    '''

    isolate_single_snps(src_snps=f"/mnt/data_3/hallep/nagnag/snps/found/{db}_nagnag_snps.txt",
                        src_sites=f"/mnt/data_3/hallep/nagnag/snps/found/{db}_nagnag_snp_containing_ssites.txt",
                        dst_snps=f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snps.txt",
                        dst_sites=f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_containing_ssites.txt",
                        fn=single_nagnag_snps)

''' ===== MINOR ALLELE FREQUENCY ===== '''

# dbSNP MINOR ALLELE FREQUENCIES
def maf_NCBI(snp:pd.Series=None) -> tuple[float, float] | list[str]:

    ''' Get minor allele frequencies from NCBI (dbSNP, ClinVar)
    
    if {snp} is not specified, get column labels:
    * "MAF_TGP"
    * "MAF_TOPMED"

    Parameters
    ----------
    snp : pandas.Series (default = None)
        dbSNP information for a particular SNP \\
        if None, return column labels

    Returns
    -------
    caf : float
        minor allele frequency based on 1000Genomes
        if not specified, -1
    topmed : float
        minor allele frequency based on TOPMED \\
        if not specified, -1
    _ : list of str
        column labels \\
        return if {snp} == None
    '''
    
    # column labels
    if type(snp) == type(None):
        return ["maf_TGP", "maf_TOPMed"]

    # alt allele index
    a = snp["alt_ind"] if "alt_ind" in snp.index else 0

    # CAF    
    if snp["CAF"] == ".":
        caf = -1
    else:
        caf = snp["CAF"].split(",")[a+1]
        caf = -1 if caf == "." else float(caf)
    
    # TOPMED
    if snp["TOPMED"] == ".":
        topmed = -1
    else:
        topmed = snp["TOPMED"].split(",")[a+1]
        topmed = -1 if topmed == "." else float(topmed)

    return caf, topmed

# GET MINOR ALLELE FREQUENCIES
def get_maf(snp_txt:str, maf_func):

    ''' Add minor allele frequency
    
    Parameters
    ----------
    snp_txt : str
        filepath to .txt file of parsed dbSNPs
    maf_func : function
        takes a row and returns that SNP's minor allele frequency/frequencies \\
        if no row passed, returns a list of column headers
    '''

    # load SNPs
    snps = pd.read_csv(snp_txt, sep="\t", index_col=0, dtype=dtypes)

    # minor allele frequencies
    frequencies = zip(*progress_map(maf_func, [s[1] for s in snps.iterrows()], n_cpu=24))

    # column labels
    columns = maf_func()

    # add columns to DataFrame
    for col,freq in zip(columns, frequencies):
        snps[col] = freq
    snps.to_csv(snp_txt, sep="\t", index=True)

    print(snps)

''' Count '''

maf_flags = {
    "dbSNP" : ["COMMON", "G5", "G5A"],
    "ClinVar" : []
}

maf_pops = {

    "dbSNP" : {
        "TGP" : "maf_TGP",
        "TOPMed" : "maf_TOPMed"
    },

    "ClinVar" : {
        "GO-ESP" : "AF_ESP",
        "ExAC" : "AF_EXAC",
        "TGP" : "AF_TGP"
    }
}

# COUNT COMMON SNPs
def count_maf(db:str) -> pd.DataFrame:

    ''' Count the number of SNPs that are "common" using different thresholds
    
    source:
    * all 3' splice sites: /mnt/data_3/hallep/nagnag/snps/{db}_3ss_snps.txt
    * nagnags: "/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snps.txt"

    destination: /mnt/data_3/hallep/nagnag/snps/statistics/{db}_maf_frequencies.txt
        * rows: (i, j)
            * i: "total", "flag", {maf_pops[db]}
            * j: "total", {maf_flags[db]}, "MAF > x%" (for x in {freqs}), "no MAF"
        * columns: (snp, freq)
            * snp: "All 3'ss", "NAGNAGs", "NAGNAGs: Create", "NAGNAGs: Alter", "NAGNAGs: Destroy"
            * freq: "count", "prop"
        * frequencies tested: 0%, 0.01%, 0.1%, 1%, 2%, 3%, 4%, 5%, 10%, 25%, 50%, 75%, 90%, 95%, 100%

    Parameters
    ----------
    ----------
    db : str {"dbSNP", "ClinVar"}
        SNP database

    Returns
    -------
    df : pandas.DataFrame
        frequencies of SNPs for given thresholds
        * columns: "All 3'ss" ("count" and "prop"), "NAGNAGs" ("count" and "prop"), "ratio"
        * rows: "total", {flags}, "MAF > x%", "no MAF"
            * "MAF > x%", "no MAF" for both "TGP" and "TOPMed"
    '''

    # load SNPs
    ss3 = pd.read_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_3ss_snps.txt", sep="\t", index_col=0, dtype=dtypes)
    nn = pd.read_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snps.txt", sep="\t", index_col=0, dtype=dtypes)

    snps = {
        "All 3'ss" : ss3,
        "NAGNAGs" : nn,
        "NAGNAGs: Create" : nn[nn["snp_type"] == "CREATE"],
        "NAGNAGS: Alter" : nn[nn["snp_type"] == "ALTER"],
        "NAGNAGS: Destroy" : nn[nn["snp_type"] == "DESTROY"]
    }
    
    # flags and frequencies to test
    freqs = [0.0, 0.0001, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 1.0]

    # minor allele frequencies
    maf_cols = {s : [df[p] for p in maf_pops[db].values()] for s,df in snps.items()}
    flag_cols = {s : [df[f] for f in maf_flags[db]] for s,df in snps.items()}

    # counts
    counts = {s : np.array([len(df)] + 
                           [len(f[f == True]) for f in flag_cols[s]] + 
                           [c for m in maf_cols[s] for c in [len(m[m >= f]) for f in freqs] + [len(m[m == -1])]])
              for s,df in snps.items()}

    # row names
    d = {
        ("i", "i") : ["total"] + ["flag" for _ in maf_flags[db]] + [pop for pop in maf_pops[db].keys() for _ in range(len(freqs)+1)],
        ("j", "j") : ["total"] + maf_flags[db] + [n for _ in range(len(maf_pops[db])) for n in [f"MAF > {f*100}%" for f in freqs] + ["no MAF"]],
    }

    # SNPs
    for s,df in snps.items():
        d.update({
            (s, "count") : counts[s],
            (s, "prop") : counts[s] / len(df)
        })

    # create DataFrame
    df = pd.DataFrame(d).set_index(keys=[("i", "i"), ("j", "j")])
    df.to_csv(f"/mnt/data_3/hallep/nagnag/snps/statistics/{db}_maf_frequencies.txt", sep="\t", index=True)

    print("\n===== dbSNP Frequencies =====")
    print(df)

    return df

''' ===== COMMON vs. RARE SNPs ===== '''

# SPLIT dbSNP SNPs
def split_dbSNP(site_type:str):

    ''' Split dbSNP SNPs and sites into common vs. rare using the "COMMON" flag
    
    source:
    * snps: /mnt/data_3/hallep/nagnag/snps/dbSNP_{site_type}_snps.txt
    * sites: /mnt/data_3/hallep/nagnag/snps/dbSNP_{site_type}_snp_containing_ssites.txt

    destination:
    * snps:
        * common: /mnt/data_3/hallep/nagnag/snps/dbSNP_common_{site_type}_snps.txt
        * rare: /mnt/data_3/hallep/nagnag/snps/dbSNP_rare_{site_type}_snps.txt
    * sites:
        * common: /mnt/data_3/hallep/nagnag/snps/dbSNP_common_{site_type}_snp_containing_ssites.txt
        * rare: /mnt/data_3/hallep/nagnag/snps/dbSNP_rare_{site_type}_snp_containing_ssites.txt

    Parameters
    ----------
    site_type : str {"3ss", "nagnag"}
        site type of SNPs to split
    '''
    
    # SNP id column
    id_col = {
        "3ss" : "snp_id_dbSNP",
        "nagnag" : "snp_ids"
    }[site_type]

    # load SNPs and sites
    snps = pd.read_csv(f"/mnt/data_3/hallep/nagnag/snps/dbSNP_{site_type}_snps.txt", sep="\t", index_col=0, dtype=dtypes)
    sites = pd.read_csv(f"/mnt/data_3/hallep/nagnag/snps/dbSNP_{site_type}_snp_containing_ssites.txt", sep="\t", index_col=0, dtype=dtypes)

    # split SNPs
    common_snps = snps[snps["COMMON"] == True]
    rare_snps = snps[snps["COMMON"] != True]

    common = set(common_snps.index.values)
    rare = set(rare_snps.index.values)

    # split sites
    common_sites = sites[sites[id_col].str.split(",").apply(lambda inds: any(i in common for i in inds))]
    rare_sites = sites[sites[id_col].str.split(",").apply(lambda inds: any(i in rare for i in inds))]

    # save
    common_snps.to_csv(f"/mnt/data_3/hallep/nagnag/snps/dbSNP_common_{site_type}_snps.txt", sep="\t", index=True)
    rare_snps.to_csv(f"/mnt/data_3/hallep/nagnag/snps/dbSNP_rare_{site_type}_snps.txt", sep="\t", index=True)

    common_sites.to_csv(f"/mnt/data_3/hallep/nagnag/snps/dbSNP_common_{site_type}_snp_containing_ssites.txt", sep="\t", index=True)
    rare_sites.to_csv(f"/mnt/data_3/hallep/nagnag/snps/dbSNP_rare_{site_type}_snp_containing_ssites.txt", sep="\t", index=True)

    print("\n===== Common SNPs =====")
    print(common_snps)

    print("\n===== Rare SNPs =====")
    print(rare_snps)

''' ===== Statistics ===== '''

# COUNT COMMON dbSNP SNPs
def count_common_dbSNP() -> tuple[int, int, int]:

    ''' Count the number of common vs. rare dbSNP SNPs
    
    source: /mnt/data_3/hallep/reference/snp_data/dbSNP/chrom_split/dbSNP_chr#.txt

    Returns
    ------- 
    total, common, rare : int
        number of all, common, and rare SNPs (respectively)
    '''
    
    common = 0
    rare = 0

    # for each chromosome
    for c in CHROMS[:-1]:

        # load snps
        snps = pd.read_csv("/mnt/data_3/hallep/reference/snp_data/dbSNP/chrom_split/dbSNP_chr#.txt".replace("#", c), sep="\t", dtype=dtypes)

        # get COMMON flag
        flag = list(progress_starmap(get_info_value, [(i, "COMMON", bool) for i in snps["info"].values], n_cpu=24))

        common += flag.count(True)
        rare += flag.count(False)

    # total
    total = common + rare

    print(f"total: {total:,} SNPs")
    print(f"common: {common:,} SNPs")
    print(f"rare: {rare:,} SNPs")

    return total, common, rare

# SNP + SITE FREQUENCIES
def snp_frequencies(dst_snps:str="/mnt/data_3/hallep/nagnag/snps/statistics/snp_frequencies.txt",
                    dst_sites:str="/mnt/data_3/hallep/nagnag/snps/statistics/site_frequencies.txt"
                    ) -> tuple[pd.DataFrame, pd.DataFrame]:

    ''' Get statistics about the types of SNPs and SNP-containing sites there are
    
    Parameters
    ----------
    dst_snps : str (default = "/mnt/data_3/hallep/nagnag/snps/statistics/snp_frequencies.txt")
        filepath to destination .txt file of SNP frequencies
    dst_sites : str (default = "/mnt/data_3/hallep/nagnag/snps/statistics/site_frequencies.txt")
        filepath to destination .txt file of site frequencies

    Returns
    -------
    snps : pandas.DataFrame
        information about the site-affecting SNPs
        * rows: "dbSNP Common", "dbSNP Rare", "ClinVar", "HGMD Splice"
        * columns: "All 3'ss", "NAGNAG", "NAGNAG: Create:", "NAGNAG: Alter", "NAGNAG: Destroy"
    sites : pandas.DataFrame
        information about SNP-containing sites
        * rows: "All 3'ss", "NAGNAGs", "1-Off", "Canonical"
        * columns: "dbSNP Common", "dbSNP Rare", "ClinVar", "HGMD Splice"
    '''

    # Totals
    ss = {st : pd.read_csv(i[1], sep="\t", dtype=dtypes).__len__() for st,i in site_info.items()}
    _, dbC, dbR = count_common_dbSNP()
    cv = pd.read_csv(var_info["ClinVar"][2], sep="\t", dtype=dtypes).__len__()
    hs = pd.read_csv(var_info["HGMD_splice"][2], sep="\t", dtype=dtypes).__len__()

    # DataFrames
    snp_df = {st : [pd.read_csv(snp_file.replace("[SNP]", db).replace("[SITE]", st), sep="\t", dtype=dtypes) for db in snp_dbs] for st in site_types}
    site_df = {db : {st : pd.read_csv(site_file.replace("[SNP]", db).replace("[SITE]", st), sep="\t", dtype=dtypes)
                        for st in site_types} for db in snp_dbs}

    for db in snp_dbs:
        n = site_df[db]["nagnag"]
        s = n["ssite_type"]
        site_df[db].update({"1off" : n[(s == "1C") | (s == "2NC")],
                            "canon" : n[s == "2C"]})

    # Splice Site-Affecting SNP Frequencies
    snps = pd.DataFrame({
        "total" : [dbC, dbR, cv, hs],
        "All 3'ss" : [len(df) for df in snp_df["3ss"]],
        "NAGNAGs" : [len(df) for df in snp_df["nagnag"]],
        "NAGNAG: Create" : [(df["snp_type"] == "CREATE").sum() for df in snp_df["nagnag"]],
        "NAGNAG: Alter" : [(df["snp_type"] == "ALTER").sum() for df in snp_df["nagnag"]],
        "NAGNAG: Destroy" : [(df["snp_type"] == "DESTROY").sum() for df in snp_df["nagnag"]]
    }, index=["dbSNP Common", "dbSNP Rare", "ClinVar", "HGMD Splice"]).T
    snps.to_csv(dst_snps, sep="\t", index_label="site")
    
    # SNP-Containing Splice Site Frequencies
    sites = pd.DataFrame({
        "Total" : [ss["all"], ss["1off"] + ss["canon"], ss["1off"], ss["canon"]],
        "dbSNP Common" : [len(df) for df in site_df["dbSNP_common"].values()],
        "dbSNP Rare" : [len(df) for df in site_df["dbSNP_rare"].values()],
        "ClinVar" : [len(df) for df in site_df["ClinVar"].values()],
        "HGMD Splice" : [len(df) for df in site_df["HGMD_splice"].values()]
    }, index=["All 3'ss", "NAGNAGs", "1-Off", "Canonical"])
    sites.to_csv(dst_sites, sep="\t", index_label="database")
    
    print("\n===== SNPs =====")
    print(snps)

    print("\n===== 3' Splice Sites =====")
    print(sites)

    return snps, sites

''' ===== RUN ===== '''

def analyze_db(db:str):
    nagnag_effect(db)
    isolate_single_3ss_snps(db)
    isolate_single_nagnag_snps(db)

dtypes.update({
    "alt_ind" : "str",
    "ref_single" : "str",
    "alt_single" : "str",
})

dbs = ["dbSNP", "ClinVar", "HGMD_splice"]

[nagnag_effect(db) for db in dbs]

''' Single SNPs '''
[isolate_single_3ss_snps(db) for db in dbs]
[isolate_single_nagnag_snps(db) for db in dbs]

dtypes.update({
    "alt_ind" : "int",
    "ref_single" : "int",
    "alt_single" : "int",
})

''' Minor Allele Frequencies '''

get_maf(snp_txt="/mnt/data_3/hallep/nagnag/snps/dbSNP_3ss_snps.txt", maf_func=maf_NCBI)
get_maf(snp_txt="/mnt/data_3/hallep/nagnag/snps/dbSNP_nagnag_snps.txt", maf_func=maf_NCBI)

count_maf("dbSNP")
count_maf("ClinVar")

''' Common vs. Rare dbSNP SNPs '''

split_dbSNP("3ss")
split_dbSNP("nagnag")

''' SNP Frequencies '''

snp_frequencies()
