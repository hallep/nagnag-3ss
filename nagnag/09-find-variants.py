''' Find splice site-affecting variants from dbSNP, ClinVar, and HGMD Splice '''

from lib import argparse, os, pd
from utils import ROOT, log_script, log_fn, single_map, multi_map
from utils.seq import CHROMS, N
from utils.vnt import dtypes, txt_cols, all_vnt_dbs, SS, DB, EFF, filter_effect

# Variant Databases
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--num-threads", type=int, help="maximum number of parallel threads to use")
parser.add_argument("-H", "--ignore-HGMD", action="store_true", help="Do not process/analyze HGMD Splice data")
args = parser.parse_args()

vnt_dbs = all_vnt_dbs
if args.ignore_HGMD:
    vnt_dbs.remove("HGMD_splice")

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
def find_variants(site_bed:str, vnt_bed:str, intersect:str, sites:pd.DataFrame, stranded:bool,
                  filter:str=None) -> tuple[list[int], list[str], list[str], list[str], list[str]]:
    
    ''' Run bedtools intersect and determine:
        * number, ids, and positions of variants in each splice site
        * affecting variants
    '''

    # run bedtools intersect
    os.system(f"bedtools intersect -a {site_bed} -b {vnt_bed} -wa -wb {"-s" if stranded else ""} > {intersect}")

    # load results
    res = pd.read_csv(intersect, sep="\t", dtype={"i":"int", "id":"str"}, index_col=0, header=None,
                      usecols=[3, 4, 5, 9, 10], names=["c1", "s1", "e1", "i", "p", "r", "c2", "s2", "e2", "id", "a", "r2"])
    
    # filter NAGNAG-affecting variants
    if filter and not res.empty:
        res = filter_effect(res, filter, stranded)
    
    if res.empty:
        return list(sites.index), [0] * len(sites), ["."] * len(sites), ["."] * len(sites), []

    # get variants that affect each site
    num, ids, pos = zip(*single_map(parse_intersect, [res.loc[i] if (i in res.index) else pd.DataFrame() for i in sites.index],
                                    n_procs=args.num_threads))

    return list(sites.index), list(num), list(ids), list(pos), list(set(res["id"]))

# Isolate affecting variants
def isolate_variants(vnt_txt:str, vnt_ids:list[str]) -> str:

    ''' Isolate nagnag-affecting variants (vnt_ids) and return as tab-separated .txt without index/headers '''
    
    # get affecting variants
    vnts = pd.read_csv(vnt_txt, sep="\t", index_col=0).loc[vnt_ids]

    # get tab-separated string
    return vnts.to_csv(None, sep="\t", index=True, header=False)

# Get affecting variants
def get_affecting_vnts(sfx:str, site_bed:str, vnt_bed:str, src_site_txt:str, dst_site_txt:str|None,
                       src_vnt_txt:str, dst_vnt_txt:str, by_chrom:bool, stranded:bool, filter:str|None,
                       indent:int):
    
    ''' Find splice site-affecting variants; isolate variants and affected sites '''

    # load sites
    sites = pd.read_csv(src_site_txt, sep="\t", index_col=0, dtype=dtypes)

    # find intersection
    log_fn("searching variants", sub=indent)
    if by_chrom:
        
        # find variants
        results = list(zip(*[find_variants(site_bed=site_bed.replace("#", c), vnt_bed=vnt_bed.replace("#", c), 
                                           intersect=f"{ROOT}/variants/bed/intersect_{c}.bed", 
                                           sites=sites[sites["chrom"] == c], stranded=stranded, filter=filter)
                             for c in CHROMS]))

        ssi = [s for chrom in results[0] for s in chrom]
        num = [n for chrom in results[1] for n in chrom]
        ids = [i for chrom in results[2] for i in chrom]
        pos = [p for chrom in results[3] for p in chrom]
        res = {s : (n, i, p) for s,n,i,p in zip(ssi, num, ids, pos)}

    else:
        
        # find variants
        results = find_variants(site_bed=site_bed, vnt_bed=vnt_bed, intersect=f"{ROOT}/variants/bed/intersect.bed",
                                sites=sites, stranded=stranded, filter=filter)
        res = {s : (n, i, p) for s,n,i,p in zip(results[0], results[1], results[2], results[3])}

    # add columns to DataFrame
    sites[f"num_vnts_{sfx}"] = [res[i][0] for i in sites.index]
    sites[f"vnt_ids_{sfx}"] = [res[i][1] for i in sites.index]
    sites[f"vnt_pos_{sfx}"] = [res[i][2] for i in sites.index]
    
    # save sites
    sites.to_csv(src_site_txt, sep="\t", index=True)

    if dst_site_txt:
        sites[sites[f"num_vnts_{sfx}"] > 0].to_csv(dst_site_txt, sep="\t", index=True)
    
    # save site-affecting variants
    log_fn("saving variants", sub=indent)
    if by_chrom:
        with open(dst_vnt_txt, "w") as file:

            # header
            header = "\t".join(txt_cols) + "\n"
            file.write(header)

            # rows
            rows = multi_map(isolate_variants, [(src_vnt_txt.replace("#", c), n) for c,n in zip(CHROMS, results[4])],
                             n_procs=args.num_threads)
            file.writelines(rows)

    else:
        vnts = pd.read_csv(src_vnt_txt, sep="\t", index_col=0, dtype=dtypes).loc[results[4]]
        vnts.to_csv(dst_vnt_txt, sep="\t", index=True)

# Combine NAGNAG-affecting variants
def process_nagnag_vnts(create_vnt_txt:str, alter_vnt_txt:str, destroy_vnt_txt:str, dst_site_txt:str, dst_vnt_txt:str,
                        dst_event_txt:str, dst_scen_txt:str, create_sfx:str, alter_sfx:str, destroy_sfx:str, stranded:bool):

    ''' Combine all NAGNAG-affecting variants and affected sites into one file,
    split sites into variant events, split events into variant splice scenarios

    **Source:**
    * `sites/1off_splice_sites.txt`
    * `sites/nagnag_3ss.txt`
    * `sites/3cs.txt`

    `{dst_site_txt}`
    ----------------
    columns: "index", "chrom", "strand", "motif_start", "motif_end", 
    "motif_seq", "motif_iflank_3ss", "motif_eflank_3ss", "event_ind", "nsite_ind", "ssite_ind",
    "ssite_type", "num_csites", "csite_starts", "csite_pos", "csite_inds",
    "num_vnts", "vnt_ids", "vnt_pos", "vnt_effects", "vnt_coords", "ref", "alt", "vnt_motif"
    * **index** is new
    * **ssite_ind** is from `3ss.txt`
    * **nsite_ind** is from `1off_splice_sites.txt` or `nagnag_3ss.txt`
    * **event_ind** is from `{dst_event_txt}`
    
    `{dst_vnt_txt}`
    ------------
    * **ssite_inds** (*str*): indices of affected splice sites in `3ss.txt`
    * **nsite_inds** (*str*): indices of affected potential sites in `nagnag_3ss.txt` or `1off_splice_sites.txt`
    * **vsite_inds** (*str*): indices of affected sites in `{dst_site_txt}`
    * **event_inds** (*str*): indices of variant event in `{dst_event_txt}`
    * **vnt_pos** (*str*): position of variant in affected site(s)
    * **vnt_effect** (*str*): effect of variant on affected site(s)
    * **tx_strands** (*str*): transcription strand of affected site(s)

    `{dst_event_txt}`
    -----------------
    columns: "index", "chrom", "motif_start", "motif_end", "strand",
    "splice_type", "ref_seq", "alt_seq", "down_seq",
    "vnt_id", "vnt_coord", "vnt_pos", "ref", "alt", "vnt_effect",
    "vsite_ind", "nsite_ind", "ssite_ind", "csite_pos", "csite_inds", "scen_inds"
    * **index** is new
    * **ssite_ind** is from `3ss.txt`
    * **nsite_ind** is from `1off_splice_sites.txt` or `nagnag_3ss.txt`
    * **vsite_ind** is from `{dst_site_txt}`

    `variants/found/{db}_nagnag_vnt_scenarios.txt`
    ----------------------------------------------
    columns: "index", "chrom", "motif_start", "motif_end", "strand",
    "splice_type", "site_splice_type", "up_end", "rtype", "phase", "acc_nums", "num_iso",
    "up_seq", "ref_seq", "alt_seq", "down_seq",
    "vnt_id", "vnt_coord", "vnt_pos", "ref", "alt", "vnt_effect",
    "event_ind", "vsite_ind", "nsite_ind", "ssite_ind", "csite_pos", "csite_inds"
    '''
    
    # load variant-containing splice sites
    df1 = pd.read_csv(f"{ROOT}/sites/1off_splice_sites.txt", sep="\t", dtype=dtypes)
    df1 = df1[df1[f"num_vnts_{create_sfx}"] > 0]

    df2 = pd.read_csv(f"{ROOT}/sites/nagnag_3ss.txt", sep="\t", dtype=dtypes)
    df2 = df2[(df2[f"num_vnts_{alter_sfx}"] > 0) | (df2[f"num_vnts_{destroy_sfx}"] > 0)]
    
    # number of variants
    numC = df1[f"num_vnts_{create_sfx}"]
    numA = df2[f"num_vnts_{alter_sfx}"]
    numD = df2[f"num_vnts_{destroy_sfx}"]

    # load site-affecting variants
    c = pd.read_csv(create_vnt_txt, sep="\t", dtype=dtypes)
    a = pd.read_csv(alter_vnt_txt, sep="\t", dtype=dtypes)
    d = pd.read_csv(destroy_vnt_txt, sep="\t", dtype=dtypes)
    vnts = pd.concat([c, a, d]).drop_duplicates(subset="id").set_index("id")

    # variant instances
    columns = {
        "chrom" : ("chrom", "chrom"),
        "strand" : ("strand", "strand"),
        "motif_start" : ("motif_start", "ssite_start"),
        "motif_end" : ("motif_end", "ssite_end"),
        "motif_seq" : ("motif_seq", "ssite_seq"),
        "motif_iflank_3ss" : ("motif_iflank_3ss", "ssite_iflank_3ss"),
        "motif_eflank_3ss" : ("motif_eflank_3ss", "ssite_eflank_3ss"),
        "nsite_ind" : ("index", "index"),
        "ssite_ind" : ("ssite_ind", "index"),
        "ssite_type" : ("ssite_type", "ssite_type"),
        "num_csites" : ("num_csites", "num_csites"),
        "csite_starts" : ("csite_starts", "csite_starts"),
        "csite_pos" : ("csite_pos", "csite_pos"),
        "csite_inds" : ("csite_inds", "csite_inds"),
    }
    d = {c : [v for n,v in zip(numC, df1[c1]) for _ in range(n)] + 
                [v for n,v in zip(numA+numD, df2[c2]) for _ in range(n)] for c,(c1,c2) in columns.items()}
    d.update({
        "vnt_id" : [i for ids in df1[f"vnt_ids_{create_sfx}"] for i in ids.split(",")] +
                        [i for idA,idD in zip(df2[f"vnt_ids_{alter_sfx}"], df2[f"vnt_ids_{destroy_sfx}"]) for i in idA.split(",") + idD.split(",") if i != "."],
        "vnt_pos" : [p for pos in df1[f"vnt_pos_{create_sfx}"] for p in pos.split(",")] +
                        [p for pA,pD in zip(df2[f"vnt_pos_{alter_sfx}"], df2[f"vnt_pos_{destroy_sfx}"]) for p in pA.split(",") + pD.split(",") if p != "."],
        "vnt_effect" : ["CREATE" for n in numC for _ in range(n)] +
                        [e for na,nd in zip(numA, numD) for e in ["ALTER"] * na + ["DESTROY"] * nd]
    })
    insts = pd.DataFrame(d)
    insts["vnt_coord"] = insts["vnt_id"].apply(lambda x: vnts.loc[x]["pos"] - 1)

    # variants with 1 effect
    vars = insts.groupby(by="vnt_id").aggregate({"vnt_effect" : lambda x: ",".join(x)})
    vars = vars[vars["vnt_effect"].apply(lambda x: len(set(x.split(","))) == 1)]
    insts = insts[insts["vnt_id"].isin(vars.index)]
    vnts = vnts.loc[vars.index]

    # alleles
    def alt_ind(inst:pd.Series) -> int:
        pos = int(inst["vnt_pos"])
        alts = str(inst["alts"]).split(",")

        # create
        if inst["vnt_effect"] == "CREATE":
            
            # desired alt allele
            if pos in [1, 4]:
                a = "A" if (inst["strand"] == "+") or stranded else "T"
            elif pos in [2, 5]:
                a = "G" if (inst["strand"] == "+") or stranded else "C"

            return alts.index(a)

        # alter or destroy
        i = 0

        # get first single-base allele
        while (alts[i].upper() not in N) or (len(alts[i]) != 1):
            i += 1

        return i
    
    def alt_allele(inst:pd.Series) -> str:

        alts = str(inst["alts"]).split(",")
        ind = int(inst["alt_ind"])
        alt = alts[ind]

        if (inst["strand"] == "+") or stranded:
            return alt
        
        return {"A":"T", "C":"G", "G":"C", "T":"A"}[alt]

    insts["alts"] = insts["vnt_id"].apply(lambda x: vnts.loc[x]["alt"])
    insts["alt_ind"] = insts.apply(lambda x: alt_ind(x), axis=1)
    insts["ref"] = insts.apply(lambda x: x["motif_seq"][int(x["vnt_pos"])], axis=1)
    insts["alt"] = insts.apply(lambda x: alt_allele(x), axis=1)
    insts["vnt_motif"] = insts.apply(lambda x: x["motif_seq"][:int(x["vnt_pos"])] + x["alt"] + x["motif_seq"][int(x["vnt_pos"])+1:], axis=1)
    
    insts.sort_values(by=["chrom", "motif_start", "strand", "vnt_coord", "vnt_effect"], ignore_index=True, inplace=True)
    insts.insert(loc=7, column="event_ind", value=insts.index)

    # --- Variants-Affected Sites --- #
    agg_fn = {c : lambda x: list(x)[0] for c in list(insts.columns)[:15]}
    agg_fn.update({c : lambda x: ",".join(map(str, x)) for c in list(insts.columns)[-9:]})
    agg_fn["vnt_effect"] = lambda x: list(x)[0]

    sites = insts.groupby(by=["ssite_ind", "nsite_ind"]).aggregate(agg_fn)
    sites.insert(loc=15, column="num_vnts", value=sites["vnt_id"].apply(lambda x: len(x.split(","))))
    sites.rename(columns={"event_ind":"event_inds", "vnt_id":"vnt_ids", "vnt_coord":"vnt_coords",
                          "ref":"refs", "alt":"alts", "vnt_motif":"vnt_motifs"}, inplace=True)
    sites.sort_values(by=["chrom", "motif_start", "strand", "vnt_effect"], inplace=True)
    sites.to_csv(dst_site_txt, sep="\t", index_label="index")

    ind_dict = sites["event_inds"].to_dict()
    insts["vsite_ind"] = insts.apply(lambda x: ind_dict[(x["ssite_ind"], x["nsite_ind"])], axis=1)

    # --- Site-Affecting Variants --- #
    columns = {
        "event_ind" : "event_inds",
        "vsite_ind" : "vsite_inds",
        "nsite_ind" : "nsite_inds",
        "ssite_ind" : "ssite_inds",
        "vnt_pos" : "vnt_pos",
        "vnt_effect" : "vnt_effect",
        "strand" : "tx_strands",
    }
    agg_fn = {c : lambda x: ",".join(map(str, x)) for c in columns}
    vars = insts.groupby(by="vnt_id").aggregate(agg_fn)

    for k,v in columns.items():
        vnts[v] = vnts.index.to_series().apply(lambda x: vars.loc[x][k])
    vnts.sort_values(by=["chrom", "pos", "tx_strands", "vnt_effect"], inplace=True)
    vnts.to_csv(dst_vnt_txt, sep="\t", index=True)

    # --- Variant Events --- #
    events = pd.DataFrame({
        "chrom" : insts["chrom"],
        "motif_start" : insts["motif_start"],
        "motif_end" : insts["motif_end"],
        "strand" : insts["strand"],
        "splice_type" : insts["csite_pos"].apply(lambda x: {"0":"PS", "1":"DS", "0,1":"AS"}[x]),
        "ref_seq" : insts["motif_seq"],
        "alt_seq" : insts["vnt_motif"],
        "down_seq" : insts["motif_eflank_3ss"].str.slice(start=-3),
        "vnt_id" : insts["vnt_id"],
        "vnt_coord" : insts["vnt_coord"],
        "vnt_pos" : insts["vnt_pos"],
        "ref" : insts["ref"],
        "alt" : insts["alt"],
        "vnt_effect" : insts["vnt_effect"],
        "vsite_ind" : insts["vsite_ind"],
        "nsite_ind" : insts["nsite_ind"],
        "ssite_ind" : insts["ssite_ind"],
        "csite_pos" : insts["csite_pos"],
        "csite_inds" : insts["csite_inds"],
    }, index=insts.index)
    
    # --- Variant Splice Scenarios --- #
    scen_inds = {i : [] for i in events.index}

    # load cleavage sites
    csites = pd.read_csv(f"{ROOT}/sites/3cs.txt", sep="\t", index_col=0)

    # columns
    d = {
        "chrom" : [],
        "motif_start" : [],
        "motif_end" : [],
        "strand" : [],

        "splice_type" : [],
        "site_splice_type" : [],

        "up_end" : [],
        "rtype" : [],
        "phase" : [],
        "acc_nums" : [],
        "num_iso" : [],

        "up_seq" : [],
        "ref_seq" : [],
        "alt_seq" : [],
        "down_seq" : [],

        "vnt_id" : [],
        "vnt_coord" : [],
        "vnt_pos" : [],
        "ref" : [],
        "alt" : [],
        "vnt_effect" : [],
        
        "event_ind" : [],
        "vsite_ind" : [],
        "nsite_ind" : [],
        "ssite_ind" : [],
        "csite_inds" : [],
        "csite_pos" : []
    }

    x = 0

    # for each event:
    for ei, event in events.iterrows():

        # cleavage sites
        cs = {int(p) : int(c) for p,c in zip(event["csite_pos"].split(","), event["csite_inds"].split(","))}

        # splice instances: (up_end, rtype, phase, up_seq) : (acc_nums, n_iso)
        insts = {0 : [], 1 : []}
        insts.update({cp : {(u,r,p,s) : (a,n) for u,r,p,s,a,n in zip(csites.loc[ci]["up_end"].split(","),
                                                                       csites.loc[ci]["rtype"].split(","),
                                                                       csites.loc[ci]["phase"].split(","),
                                                                       csites.loc[ci]["scen_eflank_5ss"].split(","),
                                                                       csites.loc[ci]["accession"].split(","),
                                                                       csites.loc[ci]["n_isoforms"].split(","))}
                            for cp,ci in cs.items()})

        # for all splice scenarios:
        for scen in set(list(insts[0]) + list(insts[1])):
            
            # location information
            for c in ["chrom", "motif_start", "motif_end", "strand"]:
                d[c].append(event[c])
            
            # splice type
            d["site_splice_type"].append(event["splice_type"])

            # alternatively-spliced
            if (scen in insts[0]) and (scen in insts[1]):
                d["splice_type"].append("AS")
                pos = [0, 1]
            
            # proximally-spliced
            elif (scen in insts[0]):
                d["splice_type"].append("PS")
                pos = [0]

            # distally-spliced
            elif (scen in insts[1]):
                d["splice_type"].append("DS")
                pos = [1]

            # scenario information
            d["up_end"].append(scen[0])
            d["rtype"].append(scen[1])
            d["phase"].append(scen[2])
            d["acc_nums"].append(",".join([insts[p][scen][0] for p in pos]))
            d["num_iso"].append(",".join([insts[p][scen][1] for p in pos]))

            # sequences
            d["up_seq"].append(scen[3][-3:])

            for c in ["ref_seq", "alt_seq", "down_seq"]:
                d[c].append(event[c])

            # variant information
            for c in ["vnt_id", "vnt_coord", "vnt_pos", "ref", "alt", "vnt_effect"]:
                d[c].append(event[c])

            # indices
            d["event_ind"].append(ei)
            d["csite_inds"].append(",".join(map(str, [cs[p] for p in pos])))
            d["csite_pos"].append(",".join(map(str, pos)))

            for c in ["vsite_ind", "nsite_ind", "ssite_ind"]:
                d[c].append(event[c])

            scen_inds[ei].append(x)
            x += 1

    # DataFrame
    scens = pd.DataFrame(d)
    scens.to_csv(dst_scen_txt, sep="\t", index_label="index")

    # add splice scenario indices to events
    events["scen_inds"] = [",".join(map(str, inds)) for inds in scen_inds.values()]
    events.to_csv(dst_event_txt, sep="\t", index=True)

# --- Generalized Functions --- #

# 3' Splice Sites
def find_3ss_vnts(db:str):

    ''' Find 3'splice site-affecting variants (call get_affecting_vnts())'''
    
    log_fn(f"Finding 3' splice site-affecting {db.replace("_", " ")} variants")

    get_affecting_vnts(sfx=db, site_bed=SS["3ss"].bed[DB[db].by_chrom], vnt_bed=DB[db].bed,
                       src_site_txt=SS["3ss"].txt, dst_site_txt=f"{ROOT}/variants/affecting/{db}_3ss_vnt_containing_ssites.txt",
                       src_vnt_txt=DB[db].txt, dst_vnt_txt=f"{ROOT}/variants/affecting/{db}_3ss_vnts.txt",
                       by_chrom=DB[db].by_chrom, stranded=DB[db].stranded, filter=None, indent=1)

# NAGNAGs
def find_nagnag_vnts(db:str):

    ''' Find NAGNAG-affecting variants; call get_affecting_vnts() and process_nagnag_vnts() '''

    log_fn(f"Finding NAGNAG-affecting {db.replace("_", " ")} variants")
    
    def by_effect(e:str):
        
        log_fn(f"NAGNAG-{e.removesuffix("e")}ing", sub=1)

        get_affecting_vnts(sfx=f"{e}_{db}", site_bed=SS[EFF[e].site_type].bed[DB[db].by_chrom],
                           vnt_bed=DB[db].bed, src_site_txt=SS[EFF[e].site_type].txt, dst_site_txt=None,
                           src_vnt_txt=DB[db].txt, dst_vnt_txt=f"{ROOT}/variants/affecting/{db}_nagnag_{EFF[e].adj}_vnts.txt",
                           by_chrom=DB[db].by_chrom, stranded=DB[db].stranded, filter=e, indent=2)

    # find
    by_effect("create")
    by_effect("alter")
    by_effect("destroy")

    # combine
    process_nagnag_vnts(create_vnt_txt=f"{ROOT}/variants/affecting/{db}_nagnag_creating_vnts.txt",
                        alter_vnt_txt=f"{ROOT}/variants/affecting/{db}_nagnag_altering_vnts.txt",
                        destroy_vnt_txt=f"{ROOT}/variants/affecting/{db}_nagnag_destroying_vnts.txt",
                        dst_site_txt=f"{ROOT}/variants/found/{db}_nagnag_vnt_containing_ssites.txt",
                        dst_vnt_txt=f"{ROOT}/variants/found/{db}_nagnag_affecting_vnts.txt",
                        dst_event_txt=f"{ROOT}/variants/found/{db}_nagnag_vnt_events.txt",
                        dst_scen_txt=f"{ROOT}/variants/found/{db}_nagnag_vnt_scenarios.txt",
                        create_sfx=f"create_{db}", alter_sfx=f"alter_{db}", destroy_sfx=f"destroy_{db}",
                        stranded=DB[db].stranded)

# ===== RUN ===== #
if __name__ == "__main__":
    log_script("09-find-variants.py")
    
    [find_3ss_vnts(db) for db in vnt_dbs]
    [find_nagnag_vnts(db) for db in vnt_dbs]
