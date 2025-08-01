import os
from tqdm import tqdm
from parallelbar import progress_map, progress_imap, progress_starmap
import pandas as pd
from snps_init import add_info_columns, add_frequency_columns, add_allmut_columns
from snps_tools import CHROMS, dtypes, vcf_columns, info_fields, effect_info, site_info, var_info, filter_create, filter_alter, filter_destroy 

''' ===== HELPER FUNCTIONS ===== '''

# AGGREGATE BY SPECIFIC COLUMN
def aggregate(df:pd.DataFrame, col:str=None) -> pd.DataFrame:

    ''' Aggregate rows in a DataFrame based on values in a specified column

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame to aggregate
    col : str (default = None)
        column by which to aggregate rows \\
        if None, use index

    Returns
    -------
    df : pandas.DataFrame
        aggregated DataFrame
    '''

    asindex = col == None
    col = col if col else df.index

    return df.groupby(df.index, as_index=asindex).agg({c : (lambda x : ",".join(map(str, x)) if len(set(x)) > 1 else x.iloc[0]) for c in df.columns})

''' ===== FIND SNPs ===== '''

# FIND SNPS IN A SPECIFIC SITE
def snps_in_site(intersect:pd.DataFrame|pd.Series) -> tuple[int, list[int|str], list[str]]:
    
    ''' Get SNPs at specific site
    
    Parameters
    ----------
    intersect : pandas.DataFrame or pandas.Series
        results of bedtools intersect for a specific site \\
        required columns:
        * "p" : position in site
        * "id" : id of SNP
        * "o" : number of overlapping bases

    Returns
    -------
    num : int
        number of SNPs in site
    pos : list of int or str
        positions of SNPs in site \\
        if none, [-1]
    ids : list of str
        SNP ids (e.g., rs number) \\
        if none, ["."]
    '''

    # no SNPs
    if intersect.empty:
        return 0, ["."], ["."]

    # 1 position (pandas.Series)
    if isinstance(intersect, pd.Series):
        
        # intersection
        num = int(intersect["o"] > 0)
        ids = [intersect["id"]] if num else ["."]
        pos = [intersect["p"]] if num else ["."]

        return num, ids, pos

    # positions with SNPs
    num = len(intersect)
    ids = intersect["id"].to_list()
    pos = intersect["p"].to_list()

    return num, ids, pos

# FIND SNPs IN SPLICE SITES
def find_snps(site_bed:str, snp_bed:str, intersect:str, sites:pd.DataFrame,
              match_strand:bool, filter=None) -> tuple[list[int], list[int], list[list[str]],
                                                       list[list[int]], list[str], list[int]]:

    ''' Identify SNPs in splice sites

    Parameters
    ----------
    site_bed : str
        filepath to .bed file of splice sites
    snp_bed : str
        filepath to .bed file of SNPs
    intersect : str
        filepath to location to save intersection results
    sites : pandas.DataFrame
        splice sites
    match_strand : bool
        whether to force strandedness when running bedtools intersect
    filter : function (default = None)
        filters NAGNAG-affecting SNPs
        * input: pandas.DataFrame (res)
        * output: pandas.DataFrame (res)
    
    Returns
    -------
    _ : list of int
        indices from {sites}
    num : list of int
        number of SNPs in each site in {sites}
    ids : list of list of str
        list of SNP ids (e.g., rs number) for each site in {sites}
        if None, "."
    pos : list of list of int
        list of SNP positions [0-5] for each site in {sites}
    snps : list of list of str
        list of SNP ids (e.g., rs number))
    inds : list of int
        index of site in {sites} for each SNP
    '''
    
    # run bedtools intersect
    os.system(f"bedtools intersect -a {site_bed} -b {snp_bed} -wao{" -s" if match_strand else ""} > {intersect}")

    # analyze results
    res = pd.read_csv(intersect, sep="\t", dtype={"i":"str","a":"str"}, index_col=0, usecols=[3, 4, 5, 9, 10, 12],
                      names=["c1", "s1", "e1", "i", "p", "r", "c2", "s2", "e2", "id", "a", "r2", "o"])
    
    # filter NAGNAG-affecting SNPs
    if filter != None:
        res = filter(res, match_strand)
    
    # remove null intersections
    res = res[res["o"] > 0]

    num, ids, pos = zip(*progress_map(snps_in_site, [res.loc[i] if (i in res.index) else pd.DataFrame() for i in sites.index], n_cpu=24))

    # SNPs
    snps = res["id"].to_list()
    inds = res.index.to_list()

    return sites.index, num, ids, pos, snps, inds

# ISOLATE SPECIFIC SNPs
def isolate_snps(src_snp_txt:str, snp_ids:list[str], indices:list[int], id_col:str="id", info_fields:dict=None) -> str:

    ''' isolate specific SNPs and get information
    
    parse info fields if {info_fields} is specified

    Parameters
    ----------
    src_snp_txt : str
        filepath to .txt file of all SNPs
    snp_ids : list of str
        indices in {src_snp_txt} to get
    indices : list of int
        splice site indices matching snps in {snp_ids}
    id_col : str (default = "id")
        column in {src_snp_txt} with SNP id number
    info_fields : dict (default = None)
        possible fields in "info"
        * keys: field name
        * values: datatype
    parallel : bool (default = False)
        whether to run add_info_columns() in parallel

    Returns
    -------
    _ : str
        tab-separated string of SNP information with the index but without the header
    '''

    # get selected SNPs
    snps = pd.read_csv(src_snp_txt, sep="\t", index_col=id_col).loc[snp_ids]

    # parse info columns
    if info_fields:
        snps = add_info_columns(df=snps, info_fields=info_fields)

    # add site indices
    snps["nsite_ind"] = indices

    # get tab-separated string
    return snps.to_csv(None, sep="\t", index=True, header=False)

# FIND SPLICE SITE-AFFECTING SNPs
def get_affecting_snps(suffix:str, by_chrom:bool, filter, src_site_bed:str, src_snp_bed:str,
                       src_site_txt:str, src_snp_txt:str, dst_snp_txt:str, id_col:str="id",
                       match_strand:bool=False, info_fields:dict=None, dst_site_txt:str=None):

    ''' Find NAGNAG-affecting SNPs
    
    src:
    * sites: {src_site_bed}
    * snps: {src_snp_bed}

    dst:
    * {src_site_txt} \\
        added columns:
        * "num_snps_{suffix}" : int \\
            number of SNPs in each 1-off NAGNAG
        * "snp_id_{suffix}" : str \\
            comma-separated list of SNP rs numbers \\
            if None, "."
        * "snp_pos_{suffix}" : str \\
            comma-separated list of SNP positions [0-5]
    * {dst_snp_txt}
        * columns in {src_snp_txt}
        * columns in {info_fields} (if specified)
        * "nsite_ind" : int \\
            index of site in {src_site_txt}
    * {dst_site_txt}
        * all columns in {src_site_txt} (+ added ones)
        * only rows where "num_snps_{suffix}" > 0

    Parameters
    ----------
    suffix : str
        suffix to add to column names
    by_chrom : bool
        whether to analyze by chromosome
    filter : function
        function that filters found SNPs
    src_site_bed : bool
        fiepath to source splice site .bed file \\
        if by_chrom == True, # will be replaced with the chromosome name
    src_snp_bed : bool
        fiepath to source splice site .bed file \\
        if by_chrom == True, # will be replaced with the chromosome name
    src_site_txt : str
        filepath to source .txt splice site database
    src_snp_txt : str
        filepath to .txt file of all SNPs \\
        if by_chrom == True, # will be replaced with the chromosome name
    dst_snp_txt : str
        filepath to .txt file to save NAGNAG-affecting SNPs
    id_col : str (default = "id")
        column in {src_snp_txt} with SNP id number
    match_strand : bool (default = False)
        whether to force strandedness when running bedtools intersect
    info_fields : dict (default = None)
        possible fields in "info"
        * keys: field name
        * values: datatype
    dst_site_txt : str (default = None)
        filepath to destination .txt splice site database \\
        only include SNP-containing sites \\
        if None, do not save
    '''

    cols = ["id"] + vcf_columns[:2] + vcf_columns[3:]

    # load sites
    sites = pd.read_csv(src_site_txt, dtype=dtypes, sep="\t", index_col=0)    

    # find intersection
    if by_chrom:

        # analyze each chromosome separately
        results = [find_snps(site_bed=src_site_bed.replace("#", c), snp_bed=src_snp_bed.replace("#", c),
                             intersect=f"/mnt/data_3/hallep/nagnag/snps/intersect/intersect_chr{c}.bed",
                             sites=sites[sites["chrom"] == f"chr{c}"], match_strand=match_strand, filter=filter) for c in CHROMS[:-2]]

        # split results
        ssi, num, ids, pos, snps, indices = zip(*results)

        ssi = [s for c in ssi for s in c]
        num = [n for c in num for n in c]
        ids = [i for c in ids for i in c]
        pos = [p for c in pos for p in c]

        res = {s : (n, i, p) for s,n,i,p in zip(ssi, num, ids, pos)}

        # add columns to splice sites DataFrame
        sites[f"num_snps_{suffix}"] = [res[i][0] for i in sites.index]
        sites[f"snp_id_{suffix}"] = [",".join(map(str, res[i][1])) for i in sites.index]
        sites[f"snp_pos_{suffix}"] = [",".join(map(str, res[i][2])) for i in sites.index]
    
    else:

        # analyze all at once
        ssi, num, ids, pos, snps, indices = find_snps(site_bed=src_site_bed, snp_bed=src_snp_bed,
                                                      intersect="/mnt/data_3/hallep/nagnag/snps/intersect/intersect.bed",
                                                      sites=sites, match_strand=match_strand, filter=filter)

        # add columns to splice sites DataFrame
        sites[f"num_snps_{suffix}"] = num
        sites[f"snp_id_{suffix}"] = [",".join(map(str, i)) for i in ids]
        sites[f"snp_pos_{suffix}"] = [",".join(map(str, p)) for p in pos]

    # save    
    sites.to_csv(src_site_txt, sep="\t", index=True)

    print("\n===== SNP-Affected Splice Sites =====")
    print(sites[sites[f"num_snps_{suffix}"] > 0])

    # save SNP-affected sites
    if dst_site_txt:
        sites[sites[f"num_snps_{suffix}"] > 0].to_csv(dst_site_txt, sep="\t", index=True)

    # Splice Site-Affecting SNPs
    if by_chrom:

        with open(dst_snp_txt, "w") as file:

            # header
            header = "\t".join(cols + ["nsite_ind"]) + "\n"
            file.write(header)

            # rows
            rows = progress_starmap(isolate_snps, [(src_snp_txt.replace("#", c), n, i, id_col, None) for c,n,i in zip(CHROMS[:-2], snps, indices)], n_cpu=24)
            # rows = [isolate_snps(src_snp_txt=src_snp_txt.replace("#", c), snp_ids=n, indices=i,
            #                      id_col=id_col, info_fields=None) for c,n,i in zip(chroms, snps, indices)]
            file.writelines(rows)

    else:

        # load all SNPs
        all_snps = pd.read_csv(src_snp_txt, sep="\t", index_col=id_col, dtype=dtypes)

        # select SNP-affected sites
        all_snps = all_snps.loc[snps]
        
        # parse info columns
        if info_fields:
            all_snps = add_info_columns(df=all_snps, info_fields=info_fields)

        # add splice site indices + aggregate
        all_snps["nsite_ind"] = indices
        all_snps = aggregate(all_snps, col=None)

        # save
        all_snps.to_csv(dst_snp_txt, sep="\t", index=True)
    
    print("\n===== Splice Site-Affecting SNPs =====")
    print(pd.read_csv(dst_snp_txt, sep="\t"))

''' ===== PROCESS SNPs ===== '''

# COMBINE NAGNAG SNPs
def combine_nagnag_snps(src_create_snps:str, src_alter_snps:str, src_destroy_snps:str, dst_snps:str, dst_sites:str,
                        create_suffix:str, alter_suffix:str, destroy_suffix:str, info_fields:dict=None, id_col:str="id"):

    ''' Get all NAGNAG-associated SNPs
    
    src:
    * 1-off: /mnt/data_3/hallep/nagnag/snps/1off_splice_sites.txt
    * canon: /mnt/data_3/hallep/nagnag/nagnag_3ss.txt

    dst:
    * {dst_sites}: "index", "chrom", "strand", "motif_start", "motif_end", 
        "motif_seq", "motif_iflank_3ss", "motif_eflank_3ss", "nsite_ind", 
        "ssite_type", "ssite_ind", "num_csites", "csite_starts", "csite_pos", "csite_inds", 
        "num_snps", "snp_ids", "snp_pos", "snp_type"
        * "index" is new
        * "ssite_ind" is from 3ss.txt
        * "nsite_ind" is from 1off_splice_sites.txt or nagnag_3ss.txt
    * {dst_snps}: {vcf_columns} + "vsite_ind" + "ssite_ind" + {info_fields} (if specified)
        * "ssite_ind" is from 3ss.txt
        * "nsite_ind" is from nagnag_3ss.txt or 1off_splice_sites.txt
        * "vsite_ind is from {dst_sites}

    Parameters
    ----------
    src_create_snps, src_alter_snps, src_destroy_snps : str
        filepath to source .txt files of NAGNAG-creating, altering, and destroying SNPs
    dst_snps : str
        filepath to destination .txt file of SNP information (from VCF)
    dst_sites : str
        filepath to destination .txt file of SNP-containing sites
    create_suffix, alter_suffix, destroy_suffix : str
        suffix for NAGNAG-creating, altering, and destroying SNP columns
    info_fields : dict (default = None)
        possible fields in "info"
        * keys: field name
        * values: datatype
    id_col : str (default = "id")
        column in {src_create_snps}, {src_alter_snps}, and {src_destroy_snps} with SNP id number
    '''
    
    # load SNP-containing splice sites
    df1 = pd.read_csv("/mnt/data_3/hallep/nagnag/snps/1off_splice_sites.txt", sep="\t", index_col=0, dtype=dtypes)
    df1 = df1[df1[f"num_snps_{create_suffix}"] > 0]

    df2 = pd.read_csv("/mnt/data_3/hallep/nagnag/nagnag_3ss.txt", sep="\t", index_col=0, dtype=dtypes)
    df2 = df2[(df2[f"num_snps_{alter_suffix}"] > 0) | (df2[f"num_snps_{destroy_suffix}"] > 0)]

    numC = df1[f"num_snps_{create_suffix}"].to_list()
    numA = df2[f"num_snps_{alter_suffix}"].to_list()
    numD = df2[f"num_snps_{destroy_suffix}"].to_list()

    # SNP-Containing Splice Sites
    sites = pd.DataFrame({
        "chrom" : df1["chrom"].to_list() + df2["chrom"].to_list(),
        "strand" : df1["strand"].to_list() + df2["strand"].to_list(),
        "motif_start" : df1["motif_start"].to_list() + df2["ssite_start"].to_list(),
        "motif_end" : df1["motif_end"].to_list() + df2["ssite_end"].to_list(),
        "motif_seq" : df1["motif_seq"].to_list() + df2["ssite_seq"].to_list(),
        "motif_iflank_3ss" : df1["motif_iflank_3ss"].to_list() + df2["ssite_iflank_3ss"].to_list(),
        "motif_eflank_3ss" : df1["motif_eflank_3ss"].to_list() + df2["ssite_eflank_3ss"].to_list(),
        "nsite_ind" : df1.index.to_list() + df2.index.to_list(),
        "ssite_ind" : df1["ssite_ind"].to_list() + df2.index.to_list(),
        "ssite_type" : df1["ssite_type"].to_list() + df2["ssite_type"].to_list(),
        "num_csites" : df1["num_csites"].to_list() + df2["num_csites"].to_list(),
        "csite_starts" : df1["csite_starts"].to_list() + df2["csite_starts"].to_list(),
        "csite_pos" : df1["csite_pos"].to_list() + df2["csite_pos"].to_list(),
        "csite_inds" : df1["csite_inds"].to_list() + df2["csite_inds"].to_list(),
        "num_snps" : numC + [a+d for a,d in zip(numA, numD)],
        "snp_ids" : df1[f"snp_id_{create_suffix}"].to_list() +
                        [",".join([r for r in z if r != "."]) for z in zip(df2[f"snp_id_{alter_suffix}"], df2[f"snp_id_{destroy_suffix}"])],
        "snp_pos" : df1[f"snp_pos_{create_suffix}"].to_list() +
                        [",".join([p for p in z if p != "."]) for z in zip(df2[f"snp_pos_{alter_suffix}"], df2[f"snp_pos_{destroy_suffix}"])],
        "snp_types" : [",".join(["CREATE"] * n) for n in numC] +
                        [",".join(["ALTER"] * a + ["DESTROY"] * d) for a,d in zip(numA, numD)]
    })
    sites.to_csv(dst_sites, sep="\t", index_label="index")

    print("\n===== SNP-Containing Splice Sites: COMBINED =====")
    print(sites)

    # SNP Information
    # load parsed SNP .vcf information
    c = pd.read_csv(src_create_snps, sep="\t", index_col=id_col, dtype=dtypes)
    a = pd.read_csv(src_alter_snps, sep="\t", index_col=id_col, dtype=dtypes)
    d = pd.read_csv(src_destroy_snps, sep="\t", index_col=id_col, dtype=dtypes)

    # add new indices
    v1 = {str(o) : str(n) for o,n in zip(df1.index, sites.index)}
    v2 = {str(o) : str(n) for o,n in zip(df2.index, sites.iloc[len(df1):].index)}
    
    n1 = {str(o) : str(n) for o,n in zip(df1.index, df1["ssite_ind"])}
    n2 = {str(o) : str(n) for o,n in zip(df2.index, df2.index)}

    def process(df:pd.DataFrame, vmap:dict[str, str], nmap:dict[str, str]):
        df["vsite_ind"] = df["nsite_ind"].apply(lambda x: vmap[str(x)])
        df["ssite_ind"] = df["nsite_ind"].apply(lambda x: nmap[str(x)])
    
    process(c, v1, n1)
    process(a, v2, n2)
    process(d, v2, n2)

    # combine + aggregate
    snps = pd.concat([c, a, d])
    snps = aggregate(snps, col=None)

    if info_fields:
        snps = add_info_columns(df=snps, info_fields=info_fields)
        
    snps.to_csv(dst_snps, sep="\t", index=True)

    print("\n===== Splice Site-Affecting SNPs: COMBINED =====")
    print(snps)

''' GENERALIZED FUNCTIONS '''

# FIND 3' SPLICE SITE-AFFECTING SNPs
def find_3ss_snps(db:str):

    ''' Find 3'splice site-affecting SNPs
    
    call get_affecting_snps() and combine_snps()

    Parameters
    ----------
    db : str
        SNP database \\
        see keys in {snp}
    '''
    
    get_affecting_snps(suffix=db, by_chrom=var_info[db][0], filter=None,
                       src_site_bed=site_info["all"][0][var_info[db][0]],
                       src_snp_bed=f"/mnt/data_3/hallep/reference/snp_data/bed/{var_info[db][1]}",
                       src_site_txt=site_info["all"][1], src_snp_txt=var_info[db][2],
                       dst_snp_txt=f"/mnt/data_3/hallep/nagnag/snps/found/{db}_3ss_snps.txt",
                       id_col=var_info[db][3], match_strand=var_info[db][4], info_fields=var_info[db][5],
                       dst_site_txt=f"/mnt/data_3/hallep/nagnag/snps/found/{db}_3ss_snp_containing_ssites.txt")

# FIND NAGNAG-AFFECTING SNPs
def find_nagnag_snps(db:str):

    ''' Find NAGNAG-affecting SNPs
    
    call get_affecting_snps() and combine_snps()

    Parameters
    ----------
    db : str
        SNP database \\
        see keys in {snp}
    '''

    def snp_effect(e:str):

        ''' Call get_affecting_snps() for effect {e} 
        
        Parameters
        ----------
        e : str {"create", "alter", "destroy"}
            NAGNAG SNP effect (see {effect})
        '''
    
        get_affecting_snps(suffix=f"{e}_{db}", by_chrom=var_info[db][0], filter=effect_info[e][1],
                           src_site_bed=site_info[effect_info[e][0]][0][var_info[db][0]],
                           src_snp_bed=f"/mnt/data_3/hallep/reference/snp_data/bed/{var_info[db][1]}",
                           src_site_txt=site_info[effect_info[e][0]][1], src_snp_txt=var_info[db][2],
                           dst_snp_txt=f"/mnt/data_3/hallep/nagnag/snps/affecting/{db}_nagnag_{effect_info[e][2]}_snps.txt",
                           id_col=var_info[db][3], match_strand=var_info[db][4], info_fields=var_info[db][5], dst_site_txt=None)

    # find
    snp_effect("create")
    snp_effect("alter")
    snp_effect("destroy")

    # combine
    combine_nagnag_snps(src_create_snps=f"/mnt/data_3/hallep/nagnag/snps/affecting/{db}_nagnag_creating_snps.txt",
                        src_alter_snps=f"/mnt/data_3/hallep/nagnag/snps/affecting/{db}_nagnag_altering_snps.txt",
                        src_destroy_snps=f"/mnt/data_3/hallep/nagnag/snps/affecting/{db}_nagnag_destroying_snps.txt",
                        dst_snps=f"/mnt/data_3/hallep/nagnag/snps/found/{db}_nagnag_snps.txt",
                        dst_sites=f"/mnt/data_3/hallep/nagnag/snps/found/{db}_nagnag_snp_containing_ssites.txt",
                        create_suffix=f"create_{db}", alter_suffix=f"alter_{db}", destroy_suffix=f"destroy_{db}",
                        info_fields=var_info[db][5], id_col=var_info[db][3])

''' ===== RUN ===== '''

# ''' dbSNP
find_3ss_snps("dbSNP")
add_info_columns(df="/mnt/data_3/hallep/nagnag/snps/found/dbSNP_3ss_snps.txt", info_fields=info_fields["dbSNP"])
add_frequency_columns(df="/mnt/data_3/hallep/nagnag/snps/found/dbSNP_3ss_snps.txt")

find_nagnag_snps("dbSNP")
add_info_columns(df="/mnt/data_3/hallep/nagnag/snps/found/dbSNP_nagnag_snps.txt", info_fields=info_fields["dbSNP"])
add_frequency_columns(df="/mnt/data_3/hallep/nagnag/snps/found/dbSNP_nagnag_snps.txt")
# '''

# ''' ClinVar
find_3ss_snps("ClinVar")
find_nagnag_snps("ClinVar")
# '''

# ''' HGMD
find_3ss_snps("HGMD_splice")
add_allmut_columns("/mnt/data_3/hallep/nagnag/snps/found/HGMD_splice_3ss_snps.txt")

find_nagnag_snps("HGMD_splice")
add_allmut_columns("/mnt/data_3/hallep/nagnag/snps/found/HGMD_splice_nagnag_snps.txt")
# '''
