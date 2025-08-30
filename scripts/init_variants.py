''' Identify all potential NAGNAGs and create .bed files for sites and variants '''

from utils import ROOT
from utils.lib import tqdm, re, pd
from utils.seq import CHROMS, get_ss_seq
from utils.vnt import txt_cols, vcf_cols, dtypes

# --- Splice Sites --- #

# Make .bed file for all 3' splice sites
def create_3ss_bed():

    ''' Create .bed files for all (uppercase) 3' splice sites

    **Source:** `sites/3ss.txt`

    `sites/3ss_uppercase.txt`
    -------------------------
    all canonical uppercase 3' splice sites

    `variants/bed/3ss.bed`
    ----------------------
    columns: "chrom", "start", "end", "name", "score", "strand"
    * **name** (*int*): index in `3ss.txt`
    * **score** (*str*): relative position
        * upstream: "u-3", "u-2", "u-1"
        * motif: 0, 1, 2, ...
        * downstream: "d1", "d2", "d3"

    `variants/bed/3ss_chr{#}.bed`
    -----------------------------
    same as `3ss.bed`, but split by chromosomes
    '''

    print("creating 3ss.bed...", end="", flush=True)

    # load 3' splice sites
    ss3 = pd.read_csv(f"{ROOT}/sites/3ss.txt", sep="\t", index_col=0)
    ss3 = ss3[ss3["uppercase"] == 1]

    # save uppercase sites
    ss3.to_csv(f"{ROOT}/sites/3ss_uppercase.txt", sep="\t", index=True)

    start = ss3[["ssite_start", "ssite_end"]].min(axis=1).values - 3
    end = ss3[["ssite_start", "ssite_end"]].max(axis=1).values + 3

    n_pos = (ss3["num_tri"].values * 3) + 6

    # .bed file
    bed = pd.DataFrame({
        "chrom" : [c for c,n in zip(ss3["chrom"], n_pos) for _ in range(n)],
        "start" : [c for s,e in zip(start, end) for c in range(s, e)],
        "end" : [c+1 for s,e in zip(start, end) for c in range(s, e)],
        "name" : [i for i,n in zip(ss3.index, n_pos) for _ in range(n)],
        "score" : [p for n in n_pos for p in ["u-3", "u-2", "u-1"] + list(range(n-6)) + ["d1", "d2", "d3"]],
        "strand" : [s for s,n in zip(ss3["strand"], n_pos) for _ in range(n)]
    })
    bed.to_csv(f"{ROOT}/variants/bed/3ss.bed", sep="\t", index=False, header=False)
    
    # split by chromosome
    for c in CHROMS:
        bed[bed["chrom"] == c].to_csv(f"{ROOT}/variants/bed/3ss_{c}.bed", sep="\t", index=False, header=False)

    print("done")

# Find 1-off NAGNAGs
def identify_1off_nagnags():

    ''' Identify 1-off NAGNAGs

    **Source:** `sites/3ss.txt`

    `sites/1off_splice_sites.txt`
    -----------------------------
    columns:  "chrom", "strand", "motif_start", "motif_end", "motif_seq", \\
        "motif_iflank_3ss", "motif_eflank_3ss", "ssite_type", "ssite_ind", \\
        "num_csites", "csite_starts", "csite_pos", "csite_inds", "off_pos", "off_coord"
    * **ssite_ind** (*int*): index of the corresponding splice site in `3ss.txt`
    * **off_pos** (*int*): position [0-5] of the "off" base within the motif
    * **off_coord** (*int*): 0-index genomic coordinate of the "off" base
    '''

    print("identifying 1-off NAGNAGs...", end="", flush=True)

    # load splice sites
    ssites = pd.read_csv(f"{ROOT}/sites/3ss.txt", sep="\t", index_col=0)

    # Canonical 1-NAGs
    c1 = ssites[(ssites["ssite_type"] == "1C") & (ssites["uppercase"] == 1)]

    inds1 = []
    starts = []
    ends = []
    csites = []
    cspos = []
    pos = []
    coord = []
    mseq = []
    iflank = []
    eflank = []

    # for each splice site:
    for i,c,cs,strand,sseq,if3,ef3 in zip(c1.index, c1["chrom"], c1["csite_starts"].astype(int), c1["strand"], c1["ssite_seq"], c1["ssite_iflank_3ss"], c1["ssite_eflank_3ss"]):
        
        xn = "((.[^A]G)|(.A[^G])|(.[^A][^G]))"
        nbg = ".[^A]G"
        nah = ".A[^G]"
        
        # NBGNAG
        if re.fullmatch(f"{xn}{nbg}", if3[-6:], re.IGNORECASE):
            inds1.append(i)       
            csites.append(cs)
            cspos.append(1)
            pos.append(1)
            mseq.append(if3[-3:] + sseq)

            if strand == "+":
                starts.append(cs-6)
                ends.append(cs)
                coord.append(cs-5)
            elif strand == "-":
                starts.append(cs+6)
                ends.append(cs)
                coord.append(cs+4)

            iflank.append(get_ss_seq(c, starts[-1], strand, 100, 0))
            eflank.append(get_ss_seq(c, ends[-1], strand, 0, 100))

        # NAHNAG
        if re.fullmatch(f"{xn}{nah}", if3[-6:], re.IGNORECASE):
            inds1.append(i)
            csites.append(cs)
            cspos.append(1)
            pos.append(2)
            mseq.append(if3[-3:] + sseq)
        
            if strand == "+":
                starts.append(cs-6)
                ends.append(cs)
                coord.append(cs-4)
            elif strand == "-":
                starts.append(cs+6)
                ends.append(cs)
                coord.append(cs+3)

            iflank.append(get_ss_seq(c, starts[-1], strand, 100, 0))
            eflank.append(get_ss_seq(c, ends[-1], strand, 0, 100))

        # NAGNBG
        if re.fullmatch(f"{nbg}{xn}", ef3[:6], re.IGNORECASE):
            inds1.append(i)
            csites.append(cs)
            cspos.append(0)
            pos.append(4)
            mseq.append(sseq + ef3[:3])
        
            if strand == "+":
                starts.append(cs-3)
                ends.append(cs+3)
                coord.append(cs+1)
            elif strand == "-":
                starts.append(cs+3)
                ends.append(cs-3)
                coord.append(cs-2)

            iflank.append(get_ss_seq(c, starts[-1], strand, 100, 0))
            eflank.append(get_ss_seq(c, ends[-1], strand, 0, 100))

        # NAGNAH
        if re.fullmatch(f"{nah}{xn}", ef3[:6], re.IGNORECASE):
            inds1.append(i)
            csites.append(cs)
            cspos.append(0)
            pos.append(5)
            mseq.append(sseq + ef3[:3])

            if strand == "+":
                starts.append(cs-3)
                ends.append(cs+3)
                coord.append(cs+2)
            elif strand == "-":
                starts.append(cs+3)
                ends.append(cs-3)
                coord.append(cs-3)

            iflank.append(get_ss_seq(c, starts[-1], strand, 100, 0))
            eflank.append(get_ss_seq(c, ends[-1], strand, 0, 100))

    # Non-Canonical 2-NAGs
    n2 = ssites[(ssites["ssite_type"] == "2NC") & (ssites["uppercase"] == 1)]

    inds2 = []

    # for each splice site:
    for i,start,strand,sseq in zip(n2.index, n2["ssite_start"], n2["strand"], n2["ssite_seq"]):

        # NBGNAG
        if re.fullmatch(".[^A]G.AG", sseq, re.IGNORECASE):
            inds2.append(i)
            pos.append(1)
            
            if strand == "+":
                coord.append(start+1)
            elif strand == "-":
                coord.append(start-2)

        # NAHNAG
        if re.fullmatch(".A[^G].AG", sseq, re.IGNORECASE):
            inds2.append(i)
            pos.append(2)

            if strand == "+":
                coord.append(start+2)
            elif strand == "-":
                coord.append(start-3)

        # NAGNBG
        if re.fullmatch(".AG.[^A]G", sseq, re.IGNORECASE):
            inds2.append(i)
            pos.append(4)

            if strand == "+":
                coord.append(start+4)
            elif strand == "-":
                coord.append(start-5)

        # NAGNAH
        if re.fullmatch(".AG.A[^G]", sseq, re.IGNORECASE):
            inds2.append(i)
            pos.append(5)

            if strand == "+":
                coord.append(start+5)
            elif strand == "-":
                coord.append(start-6)

    # DataFrame
    df = pd.DataFrame({
        "chrom" : [c1["chrom"][i] for i in inds1] + [n2["chrom"][i] for i in inds2],
        "strand" : [c1["strand"][i] for i in inds1] + [n2["strand"][i] for i in inds2],
        "motif_start" : starts + [n2["ssite_start"][i] for i in inds2],
        "motif_end" : ends + [n2["ssite_end"][i] for i in inds2],
        "motif_seq" : mseq + [n2["ssite_seq"][i] for i in inds2],
        "motif_iflank_3ss" : iflank + [n2["ssite_iflank_3ss"][i] for i in inds2],
        "motif_eflank_3ss" : eflank + [n2["ssite_eflank_3ss"][i] for i in inds2],
        "ssite_type" : [c1["ssite_type"][i] for i in inds1] + [n2["ssite_type"][i] for i in inds2],
        "ssite_ind" : inds1 + inds2,
        "num_csites" : [1] * len(inds1) + [n2["num_csites"][i] for i in inds2],
        "csite_starts" : csites + [n2["csite_starts"][i] for i in inds2],
        "csite_pos" : cspos + [n2["csite_pos"][i] for i in inds2],
        "csite_inds" : [c1["csite_inds"][i] for i in inds1] + [n2["csite_inds"][i] for i in inds2],
        "off_pos" : pos,
        "off_coord" : coord,
    })
    df.to_csv(f"{ROOT}/sites/1off_splice_sites.txt", sep="\t", index_label="index")

    print("done")

# Make .bed file for (potential) NAGNAGs
def create_nagnag_1off_bed():

    ''' Create .bed files for NAGNAGs

    **Source:**
    * `sites/1off_splice_sites.txt`
    * `sites/nagnag_3ss.txt`

    `variants/bed/1off.bed`
    -----------------------
    columns: "chrom", "start", "end", "name", "score", "strand"
    * **name** (*int*): index in `1off_splice_sites.txt`
    * **score** (*int*): relative position [0-5]

    `variants/bed/canon.bed`
    ------------------------
    columns: "chrom", "start", "end", "name", "score", "strand"
    * **name** (*int*): index in `nagnag_3ss.txt`
    * **score** (*int*): relative position [0-5]
    '''

    print("creating 1off.bed and canon.bed...", end="", flush=True)

    # 1-off NAGNAGs
    off = pd.read_csv(f"{ROOT}/sites/1off_splice_sites.txt", sep="\t", index_col=0, dtype=dtypes)

    # DataFrame
    bed = pd.DataFrame({
        "chrom" : off["chrom"],
        "start" : off["off_coord"],
        "end" : off["off_coord"]+1,
        "name" : off.index,
        "score" : off["off_pos"],
        "strand" : off["strand"]
    })
    bed.to_csv(f"{ROOT}/variants/bed/1off.bed", sep="\t", index=False, header=False)

    for c in CHROMS:
        bed[bed["chrom"] == c].to_csv(f"{ROOT}/variants/bed/1off_{c}.bed", sep="\t", index=False, header=False)

    # Canonical NAGNAGs
    canon = pd.read_csv(f"{ROOT}/sites/nagnag_3ss.txt", sep="\t", index_col=0)

    # DataFrame
    bed = pd.DataFrame({
        "chrom" : [c for c in canon["chrom"] for _ in range(6)],
        "start" : [s+p if (r == "+") else s-p-1 for s,r in zip(canon["ssite_start"],canon["strand"]) for p in range(6)],
        "end" : [s+p+1 if (r == "+") else s-p for s,r in zip(canon["ssite_start"],canon["strand"]) for p in range(6)],
        "name" : [i for i in canon.index for _ in range(6)],
        "score" : [s for _ in canon.index for s in range(6)],
        "strand" : [s for s in canon["strand"] for _ in range(6)]
    })
    bed.to_csv(f"{ROOT}/variants/bed/canon.bed", sep="\t", index=False, header=False)

    for c in CHROMS:
        bed[bed["chrom"] == c].to_csv(f"{ROOT}/variants/bed/canon_{c}.bed", sep="\t", index=False, header=False)

    print("done")

# --- Variants --- '''

# format chromosome
def fmt_chrom(c:str) -> str:
    ''' Format variant chromosome:
    * add 'chr' prefix
    * replace 'MT' with 'M' (to match RefSeq)
    '''

    c = c.upper().removeprefix("CHR")
    
    if c == "MT":
        return "chrM"

    return f"chr{c}"

# .vcf to .txt file to .bed file
def vcf2txt2bed(vcf:str, txt:str, bed:str, by_chrom:bool=False, split_freq:bool=False):
    
    ''' Convert .vcf file to .txt file, and .txt file to .bed file
    
    `{txt}`
    -------
    columns: "id", "chrom", "pos", "ref", "alt", "qual", "filter", "info"
    * **pos** (*int*): 1-index coordinate

    `{bed}`
    -------
    columns: "chrom", "start", "end", "name", "score", "strand"
    * **name** (*str*): rs accession number (from "id" column)
    * **score** (*str*): alternate base(s) (from "alt" column)
    * **strand** (*str*): always "." (unspecified)
    '''

    if by_chrom:

        print(f"converting {vcf.split("/")[-1]} to {txt.split("/")[-1]} and {bed.split("/")[-1]}...")

        # open files
        if split_freq:
            freq = {True : "common", False : "rare"}
            txt_files = {c : {f : open(txt.replace("#", c).replace("$", l), "w") for f,l in freq.items()} for c in CHROMS}
            bed_files = {c : {f : open(bed.replace("#", c).replace("$", l), "w") for f,l in freq.items()} for c in CHROMS}

        else:
            txt_files = {c : {True : open(txt.replace("#", c), "w")} for c in CHROMS}
            bed_files = {c : {True : open(bed.replace("#", c), "w")} for c in CHROMS}

        # header
        header = "\t".join(txt_cols) + "\n"
        [f.write(header) for c in txt_files.values() for f in c.values()]

        # for each line:
        for line in tqdm(open(vcf, "r")):

            # if variant (not comment):
            if line[0] != "#":
                                    
                # parse line
                l = line.strip().split("\t")
                chrom, pos, rs, ref, alt = fmt_chrom(l[0]), int(l[1]), l[2], l[3], l[4]

                # if single-base substitution on canonical chromosome:
                if (chrom in CHROMS) and (len(ref) == 1) and (any((len(a) == 1) for a in alt.split(","))):
                    
                    com = "COMMON=1" in l[7].split(";") if split_freq else True

                    # .txt file
                    txt_files[chrom][com].write("\t".join([rs] + l[:2] + l[3:]) + "\n")

                    # .bed file
                    bed_files[chrom][com].write("\t".join([chrom, str(pos-1), str(pos), rs, alt, "."]) + "\n")

        # close files
        [f.close() for c in txt_files.values() for f in c.values()]
        [f.close() for c in bed_files.values() for f in c.values()]

    else:

        print(f"converting {vcf.split("/")[-1]} to {txt.split("/")[-1]} and {bed.split("/")[-1]}...", end="", flush=True)

        # .txt file
        df = pd.read_csv(vcf, sep="\t", comment="#", header=None, names=vcf_cols,
                         dtype={0:"str", 1:"int", 2:"str", 3:"str", 4:"str", 5:"str", 6:"str", 7:"str"})
        df = df[["id", "chrom", "pos", "ref", "alt", "qual", "filter", "info"]]

        # add "chr" to chromosome name
        df["chrom"] = df["chrom"].apply(lambda x: fmt_chrom(x))

        # filter single-base substitutions on canonical chromosomes
        df = df[(df["chrom"].apply(lambda x: x in CHROMS + ["MT"])) & (df["ref"].str.len() == 1) &
                (df["alt"].apply(lambda x: any((len(a) == 1) for a in x.split(","))))]

        # .bed file
        def txt2bed(df:pd.DataFrame, bed_file:str):
            df["start"] = df["pos"] - 1
            df["strand"] = ["."] * len(df)
            df = df[["chrom", "start", "pos", "id", "alt", "strand"]]
            df.to_csv(bed_file, sep="\t", index=False, header=False)

        if split_freq:

            # common
            dfC = df[df["info"].apply(lambda x: "COMMON=1" in x.split(";"))]
            dfC.to_csv(txt.replace("$", "common"))
            txt2bed(dfC, bed.replace("$", "common"))

            # rare
            dfR = df[df["info"].apply(lambda x: "COMMON=1" not in x.split(";"))]            
            dfR.to_csv(txt.replace("$", "rare"))
            txt2bed(dfR, bed.replace("$", "rare"))

        else:

            # all
            df.to_csv(txt, sep="\t", index=False)
            txt2bed(df, bed)

    print("done")

# SQL .txt file to .bed file
def sql2txt2bed(sql:str, txt:str, bed:str, coords:str):

    ''' Convert a .txt SQL to a .bed file and save copy of .txt file 
    
    `{txt}`
    -------
    columns: "id", "chrom", "pos", "pos_end", "strand", "ref", "alt", "gene", "disease"
    * **pos** (*int*): 1-index coordinate of mutation start
    * **pos_end** (*int*): 1-index coordinate of mutation end (inclusive)

    `{bed}`
    -------
    columns: "chrom", "start", "end", "name", "score", "strand"
    * **name** (*str*): accession number (from "acc_num" column)
    * **score** (*str*): alternate base(s) (from "base" column)
    '''

    print(f"converting {sql.split("/")[-1]} to {txt.split("/")[-1]} and {bed.split("/")[-1]}...", end="", flush=True)

    # load HGMD coordinates
    hg38 = pd.read_csv(coords, sep="\t", index_col="acc_num")
    hg38["chromosome"] = hg38["chromosome"].apply(lambda x: fmt_chrom(x))
    valid = set(hg38.index)

    # load SQL table
    table = pd.read_csv(sql, sep="\t")
    table = table[table["acc_num"].apply(lambda x: x in valid)]
    ids = table["acc_num"]

    # .txt file
    df = pd.DataFrame({
        "id" : ids,
        "chrom" : [hg38.loc[i]["chromosome"] for i in ids],
        "pos" : [hg38.loc[i]["coordSTART"] for i in ids],
        "pos_end" : [hg38.loc[i]["coordEND"] for i in ids],
        "strand" : [hg38.loc[i]["strand"] for i in ids],
        "ref" : table["base"].str.split("[->]", regex=True).apply(lambda x: x[0]),
        "alt" : table["base"].str.split("[->]", regex=True).apply(lambda x: x[1]),
        "gene" : table["gene"],
        "disease" : table["disease"]
    })

    # filter single-base substitutions on canonical chromosomes
    df = df[(df["chrom"].apply(lambda x: x in CHROMS)) & (df["ref"].str.len() == 1) &
            (df["alt"].apply(lambda x: any((len(a) == 1) for a in x.split(","))))]

    df.to_csv(txt, sep="\t", index=False)

    # .bed file
    df["pos"] = df["pos"] - 1
    df = df[["chrom", "pos", "pos_end", "id", "alt", "strand"]]
    df.to_csv(bed, sep="\t", index=False, header=False)

    print("done")

# --- Run --- #

create_3ss_bed()
identify_1off_nagnags()
create_nagnag_1off_bed()

vcf2txt2bed(vcf=f"{ROOT}/src/dbSNP.vcf", txt=f"{ROOT}/variants/src/dbSNP_$_#.txt",
            bed=f"{ROOT}/variants/bed/dbSNP_$_#.bed", by_chrom=True, split_freq=True)

vcf2txt2bed(vcf=f"{ROOT}/src/ClinVar.vcf", txt=f"{ROOT}/variants/src/ClinVar.txt",
            bed=f"{ROOT}/variants/bed/ClinVar.bed", by_chrom=False, split_freq=False)

sql2txt2bed(sql=f"{ROOT}/src/HGMD_pro/splice.txt", txt=f"{ROOT}/variants/src/HGMD_splice.txt",
            bed=f"{ROOT}/variants/bed/HGMD_splice.bed", coords=f"{ROOT}/src/HGMD_pro/hg38_coords.txt")
