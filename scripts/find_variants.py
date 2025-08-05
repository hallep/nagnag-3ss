from lib import ROOT, os, progress_map, progress_starmap, pd
from sequence import CHROMS
from scripts.variants import dtypes, txt_cols, SS, DB, EFF, filter_effect

# Parse bedtools intersect output
def parse_intersect(intersect:pd.DataFrame) -> tuple[int, str, str]:

    ''' Get variants affecting specific site '''

    # no variants
    if intersect.empty:
        return 0, ".", "."
    
    # 1 variant
    if isinstance(intersect, pd.Series):
        return 1, intersect["id"], intersect["p"]

    # multiple variants
    num = len(intersect)
    ids = ",".join(map(str, intersect["id"].to_list()))
    pos = ",".join(map(str, intersect["p"].to_list()))

    return num, ids, pos

# Find affecting variants
def find_variants(site_bed:str, vnt_bed:str, intersect:str, sites:pd.DataFrame, match_strand:bool,
                  filter:str=None) -> tuple[list[int], list[str], list[str], list[str], list[str], list[str]]:
    
    ''' Run bedtools intersect and determine:
        * number, ids, and positions of variants in each splice site
        * affecting variants and affected sites' indices
    '''

    # run bedtools intersect
    os.system(f"bedtools intersect -a {site_bed} -b {vnt_bed} -wa -wb {"-s" if match_strand else ""} > {intersect}")

    # load results
    res = pd.read_csv(intersect, sep="\t", dtype={"i":"int","a":"str"}, index_col=0, header=None, usecols=[3, 4, 5, 9, 10],
                      names=["c1", "s1", "e1", "i", "p", "r", "c2", "s2", "e2", "id", "a", "r2"])

    # filter NAGNAG-affecting variants
    if filter:
        res = filter_effect(res, filter, match_strand)
    
    # get variants that affect each site
    num, ids, pos = zip(*progress_map(parse_intersect, [res.loc[i] if (i in res.index) else pd.DataFrame() for i in sites.index], n_cpu=24))

    # get sites that affect each variant
    res["vnts"] = res.index.to_series()
    agg = res.groupby(by="id").aggregate({"vnts" : lambda x: ",".join(map(str, x))})

    vnt_ids = agg.index.to_list()
    site_inds = agg["vnts"].to_list()

    return sites.index, list(num), list(ids), list(pos), vnt_ids, site_inds

# Isolate affecting variants
def isolate_variants(vnt_txt:str, vnt_ids:list[str], site_inds:list[str]) -> str:

    ''' Isolate nagnag-affecting variants (vnt_ids), add affected sites' indices (site_inds),
        and return as tab-separated .txt without index/headers '''
    
    # get affecting variants
    vnts = pd.read_csv(vnt_txt, sep="\t", index_col=0).loc[vnt_ids]

    # add site indices
    vnts["nsite_ind"] = site_inds

    # get tab-separated string
    return vnts.to_csv(None, sep="\t", index=True, header=False)

# Get affecting variants
def get_affecting_vnts(sfx:str, site_bed:str, vnt_bed:str, src_site_txt:str, dst_site_txt:str|None,
                       src_vnt_txt:str, dst_vnt_txt:str, by_chrom:bool, match_strand:bool, filter:str|None):
    
    ''' Find splice site-affecting variants '''

    # load sites
    sites = pd.read_csv(src_site_txt, sep="\t", index_col=0, dtype=dtypes)

    # find intersection
    if by_chrom:

        # find variants
        results = zip(*[find_variants(site_bed=site_bed.replace("#", c), vnt_bed=vnt_bed.replace("#", c),
                                      intersect=f"{ROOT}/variants/bed/intersect_{c}.bed",
                                      sites=sites[sites["chrom"] == c], match_strand=match_strand, filter=filter)
                        for c in CHROMS])

        ssi = [s for chrom in results[0] for s in chrom]
        num = [n for chrom in results[1] for n in chrom]
        ids = [i for chrom in results[2] for i in chrom]
        pos = [p for chrom in results[3] for p in chrom]
        res = {s : (n, i, p) for s,n,i,p in zip(ssi, num, ids, pos)}

        # add columns to DataFrame
        sites[f"num_vnts_{sfx}"] = [res[i][0] for i in sites.index]
        sites[f"vnt_ids_{sfx}"] = [res[i][1] for i in sites.index]
        sites[f"vnt_pos_{sfx}"] = [res[i][2] for i in sites.index]

    else:
        
        # find variants
        results = find_variants(site_bed=site_bed, vnt_bed=vnt_bed, intersect=f"{ROOT}/variants/bed/intersect.bed",
                                sites=sites, match_strand=match_strand, filter=filter)

        # add columns to DataFrame
        sites[f"num_vnts_{sfx}"] = results[1]
        sites[f"vnt_ids_{sfx}"] = results[2]
        sites[f"vnt_pos_{sfx}"] = results[3]

    # save sites
    sites.to_csv(src_site_txt, sep="\t", index=True)

    if dst_site_txt:
        sites[sites[f"num_vnts_{sfx}"] > 0].to_csv(dst_site_txt, sep="\t", index=True)

    # save site-affecting variants
    if by_chrom:
        with open(dst_vnt_txt, "w") as file:

            # header
            header = "\t".join(txt_cols + ["nsite_ind"]) + "\n"
            file.write(header)

            # rows
            rows = progress_starmap(isolate_variants, [(src_vnt_txt.replace("#", c), n, i, None)
                                                       for c,n,i in zip(CHROMS, results[4], results[5])], n_cpu=24)
            file.writelines(rows)

    else:
        vnts = pd.read_csv(src_vnt_txt, sep="\t", index_col=0, dtype=dtypes).loc[results[4]]
        vnts["nsite_ind"] = results[5]
        vnts.to_csv(dst_vnt_txt, sep="\t", index=True)

# Combine NAGNAG-affecting variants
def combine_nagnag_vnts(create_vnt_txt:str, alter_vnt_txt:str, destroy_vnt_txt:str, dst_vnt_txt:str,
                        dst_site_txt:str, create_sfx:str, alter_sfx:str, destroy_sfx:str):

    # load Variant-containing splice sites
    df1 = pd.read_csv(f"{ROOT}/sites/1off_splice_sites.txt", sep="\t", index_col=0, dtype=dtypes)
    df1 = df1[df1[f"num_vnts_{create_sfx}"] > 0]

    df2 = pd.read_csv(f"{ROOT}/sites/nagnag_3ss.txt", sep="\t", index_col=0, dtype=dtypes)
    df2 = df2[(df2[f"num_vnts_{alter_sfx}"] > 0) | (df2[f"num_vnts_{destroy_sfx}"] > 0)]

    numC = df1[f"num_vnts_{create_sfx}"].to_list()
    numA = df2[f"num_vnts_{alter_sfx}"].to_list()
    numD = df2[f"num_vnts_{destroy_sfx}"].to_list()

    # Variant-Containing Splice Sites
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
        "num_vnts" : numC + [a+d for a,d in zip(numA, numD)],
        "vnt_ids" : df1[f"vnt_ids_{create_sfx}"].to_list() +
                        [",".join([r for r in z if r != "."]) for z in zip(df2[f"vnt_ids_{alter_sfx}"], df2[f"vnt_ids_{destroy_sfx}"])],
        "vnt_pos" : df1[f"vnt_pos_{create_sfx}"].to_list() +
                        [",".join([p for p in z if p != "."]) for z in zip(df2[f"vnt_pos_{alter_sfx}"], df2[f"vnt_pos_{destroy_sfx}"])],
        "vnt_effects" : [",".join(["CREATE"] * n) for n in numC] +
                        [",".join(["ALTER"] * a + ["DESTROY"] * d) for a,d in zip(numA, numD)]
    })
    sites.to_csv(dst_site_txt, sep="\t", index_label="index")

    # Variant Information
    # load parsed Variant .vcf information
    c = pd.read_csv(create_vnt_txt, sep="\t", index_col=0, dtype=dtypes)
    a = pd.read_csv(alter_vnt_txt, sep="\t", index_col=0, dtype=dtypes)
    d = pd.read_csv(destroy_vnt_txt, sep="\t", index_col=0, dtype=dtypes)

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
    vnts = pd.concat([c, a, d])

    agg_fn = {c : lambda x: ",".join(map(str, x)) if len(set(x)) > 1 else x.iloc[0] for c in vnts.columns}
    vnts = vnts.groupby(vnts.index, as_index=True).aggregate(agg_fn)
    
    vnts.to_csv(dst_vnt_txt, sep="\t", index=True)

# ===== Generalized Functions ===== #

# 3' Splice Sites
def find_3ss_vnts(db:str):

    ''' Find 3'splice site-affecting variants (call get_affecting_vnts())'''
    
    get_affecting_vnts(sfx=db, site_bed=SS["3ss"].bed[DB[db].by_chrom], vnt_bed=DB[db].bed,
                       src_site_txt=SS["3ss"].txt, dst_site_txt=f"{ROOT}/variants/found/{db}_3ss_vnt_containing_ssites.txt",
                       src_vnt_txt=DB[db].txt, dst_vnt_txt=f"{ROOT}/variants/found/{db}_3ss_vnts.txt",
                       by_chrom=DB[db].by_chrom, match_strand=DB[db].match_strand, filter=None)

# FIND NAGNAG-AFFECTING SNPs
def find_nagnag_snps(db:str):

    ''' Find NAGNAG-affecting SNPs (call get_affecting_vnts() and combine_nagnag_vnts()) '''

    def snp_effect(e:str):

        get_affecting_vnts(sfx=f"{e}_{db}", site_bed=SS[EFF[e].site_type].bed[DB[db].by_chrom],
                           vnt_bed=DB[db].bed, src_site_txt=SS[EFF[e].site_type].txt, dst_site_txt=None,
                           src_vnt_txt=DB[db].txt, dst_vnt_txt=f"{ROOT}/variants/found/{db}_nagnag_{EFF[e].adj}_vnts.txt",
                           by_chrom=DB[db].by_chrom, match_strand=DB[db].match_strand, filter=e)

    # find
    snp_effect("create")
    snp_effect("alter")
    snp_effect("destroy")

    # combine
    combine_nagnag_vnts(create_vnt_txt=f"{ROOT}/variants/found/{db}_nagnag_creating_vnts.txt",
                        alter_vnt_txt=f"{ROOT}/variants/found/{db}_nagnag_altering_vnts.txt",
                        destroy_vnt_txt=f"{ROOT}/variants/found/{db}_nagnag_destroying_vnts.txt",
                        dst_vnt_txt=f"{ROOT}/variants/found/{db}_nagnag_affecting_vnts.txt",
                        dst_site_txt=f"{ROOT}/variants/found/{db}_nagnag_vnts_containing_ssites.txt",
                        create_sfx=f"create_{db}", alter_sfx=f"alter_{db}", destroy_sfx=f"destroy_{db}")

find_3ss_vnts("dbSNP")
find_nagnag_snps("dbSNP")

find_3ss_vnts("ClinVar")
find_nagnag_snps("ClinVar")

find_3ss_vnts("HGMD_splice")
find_nagnag_snps("HGMD_splice")
