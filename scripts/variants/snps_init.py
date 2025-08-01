import re
import pandas as pd
import maskpass
from sqlalchemy import create_engine, MetaData
from tqdm import tqdm
from parallelbar import progress_map, progress_imap, progress_starmap
import numpy as np
from sequence import get_ss_seq
from snps_tools import CHROMS, dtypes, vcf_columns, info_fields, info_defaults

''' ===== 3' SPLICE SITES ===== '''

# .BED FILES: ALL SPLICE SITES
def create_3ss_bed():

    ''' Create .bed files for all (uppercase) 3' splice sites

    src: /mnt/data_3/hallep/nagnag/3ss.txt

    save uppercase sites as /mnt/data_3/hallep/nagnag/3ss_uppercase.txt

    dst:
    * all: /mnt/data_3/hallep/nagnag/snps/3ss.bed
            * columns: "chrom", "start", "end", "name", "score", "strand"
            * "name" is index in 3ss.txt
            * "score" is relative position \\
                * 3 intronic bases upstream: "u-3", "u-2", "u-1"
                * motif bases: 0, 1, 2, ... 
                * 3 exonic bases downstream: "d1", "d2", "d3"
    * by chromosome: /mnt/data_3/hallep/nagnag/snps/intersect/3ss_chr{#}.bed
    '''

    # load all 3' splice sites
    ss3 = pd.read_csv("/mnt/data_3/hallep/nagnag/3ss.txt", sep="\t", index_col=0)
    ss3 = ss3[ss3["uppercase"] == 1]

    # save uppercase sites
    ss3.to_csv("/mnt/data_3/hallep/nagnag/3ss_uppercase.txt", sep="\t", index=True)

    start = ss3[["ssite_start", "ssite_end"]].min(axis=1).values - 3
    end = ss3[["ssite_start", "ssite_end"]].max(axis=1).values + 3

    n_pos = (ss3["num_tri"].values * 3) + 6

    # create .bed
    bed = pd.DataFrame({
        "chrom" : [c for c,n in zip(ss3["chrom"], n_pos) for _ in range(n)],
        "start" : [c for s,e in zip(start, end) for c in range(s, e)],
        "end" : [c+1 for s,e in zip(start, end) for c in range(s, e)],
        "name" : [i for i,n in zip(ss3.index, n_pos) for _ in range(n)],
        "score" : [p for n in n_pos for p in ["u-3", "u-2", "u-1"] + list(range(n-6)) + ["d1", "d2", "d3"]],
        "strand" : [s for s,n in zip(ss3["strand"], n_pos) for _ in range(n)]
    })
    bed.to_csv(f"/mnt/data_3/hallep/nagnag/snps/intersect/3ss.bed", sep="\t", index=False, header=False)

    print(bed)

    # save split
    for c in CHROMS[:-2]:
        bed[bed["chrom"] == f"chr{c}"].to_csv(f"/mnt/data_3/hallep/nagnag/snps/intersect/3ss_chr{c}.bed", sep="\t", index=False, header=False)

# GET 1-OFF NAGNAGs
def get_1off_nagnags():

    ''' Get 1-off NAGNAGs

    src: /mnt/data_3/hallep/nagnag/3ss.txt

    dst: /mnt/data_3/hallep/nagnag/snps/1off_splice_sites.txt \\
    "chrom", "strand", "motif_start", "motif_end", "motif_seq", \\
        "motif_iflank_3ss", "motif_eflank_3ss", "ssite_type", "ssite_ind", \\
        "num_csites", "csite_starts", "csite_pos", "csite_inds", "off_pos", "off_coord"
            * "ssite_ind" is the index of the corresponding splice site in 3ss.txt \\
                there are some duplicates!
            * "off_pos" is the position [0-5] of the "off" base within the motif
            * "off_coord" is the genomic coordinate of the "off" base
    '''

    # load splice sites
    ssites = pd.read_csv("/mnt/data_3/hallep/nagnag/3ss.txt", sep="\t", index_col=0)

    # CANONICAL 1-NAGs
    c1 = ssites[(ssites["ssite_type"] == "1C") & (ssites["uppercase"] == 1)]
    c1.loc[:, "csite_starts"] = c1["csite_starts"].astype(int)

    inds1 = []
    starts = []
    ends = []
    csites = []
    cspos = []
    pos = []
    coord = []
    seq = []
    iflank = []
    eflank = []

    # for each splice site:
    for i,c,cs,strand,sseq,if3,ef3 in zip(c1.index, c1["chrom"], c1["csite_starts"], c1["strand"], c1["ssite_seq"], c1["ssite_iflank_3ss"], c1["ssite_eflank_3ss"]):
        
        xn = "((.[^A]G)|(.A[^G])|(.[^A][^G]))"
        nbg = ".[^A]G"
        nah = ".A[^G]"
        
        # NBGNAG
        if re.fullmatch(f"{xn}{nbg}", if3[-6:], re.IGNORECASE):
            inds1.append(i)       
            csites.append(cs)
            cspos.append(1)
            pos.append(1)
            seq.append(if3[-3:] + sseq)

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
            seq.append(if3[-3:] + sseq)
        
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
            seq.append(sseq + ef3[:3])
        
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
            seq.append(sseq + ef3[:3])

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

    # NON-CANONICAL 2-NAGs
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
        "motif_seq" : seq + [n2["ssite_seq"][i] for i in inds2],
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
    df.to_csv("/mnt/data_3/hallep/nagnag/snps/1off_splice_sites.txt", sep="\t", index_label="index")

    print(df)

# .BED FILES: (POTENTIAL) NAGNAGs
def create_nagnag_1off_bed():

    ''' Create .bed files for NAGNAGs

    src:
    * 1-off: /mnt/data_3/hallep/nagnag/snps/intersect/1off_splice_sites.txt
    * canon: /mnt/data_3/hallep/nagnag/nagnag_3ss.txt

    dst:
    * 1-off: /mnt/data_3/hallep/nagnag/snps/1off.bed
            * columns: "chrom", "start", "end", "name", "score", "strand"
            * "name" is index in 1off_splice_sites.txt
            * "score" is relative position [0-5]
    * canon: /mnt/data_3/hallep/nagnag/snps/intersect/canon.bed
            * columns: "chrom", "start", "end", "name", "score", "strand"
            * "name" is index in nagnag_3ss.txt
            * "score" is relative position [0-5]
    '''

    # 1-OFF NAGNAGs
    off = pd.read_csv("/mnt/data_3/hallep/nagnag/snps/1off_splice_sites.txt",
                      sep="\t", index_col=0, dtype=dtypes)

    # DataFrame
    bed = pd.DataFrame({
        "chrom" : off["chrom"],
        "start" : off["off_coord"],
        "end" : off["off_coord"]+1,
        "name" : off.index,
        "score" : off["off_pos"],
        "strand" : off["strand"]
    })
    bed.to_csv("/mnt/data_3/hallep/nagnag/snps/intersect/1off.bed", sep="\t", index=False, header=False)

    # save split
    for c in CHROMS[:-2]:
        bed[bed["chrom"] == f"chr{c}"].to_csv(f"/mnt/data_3/hallep/nagnag/snps/intersect/1off_chr{c}.bed", sep="\t", index=False, header=False)

    print("\n===== 1-Off NAGNAGs =====")
    print(bed)

    # CANONICAL NAGNAGs
    canon = pd.read_csv("/mnt/data_3/hallep/nagnag/nagnag_3ss.txt", sep="\t", index_col=0)

    # DataFrame
    bed = pd.DataFrame({
        "chrom" : [c for c in canon["chrom"] for _ in range(6)],
        "start" : [s+p if (r == "+") else s-p-1 for s,r in zip(canon["ssite_start"],canon["strand"]) for p in range(6)],
        "end" : [s+p+1 if (r == "+") else s-p for s,r in zip(canon["ssite_start"],canon["strand"]) for p in range(6)],
        "name" : [i for i in canon.index for _ in range(6)],
        "score" : [s for _ in canon.index for s in range(6)],
        "strand" : [s for s in canon["strand"] for _ in range(6)]
    })
    bed.to_csv("/mnt/data_3/hallep/nagnag/snps/intersect/canon.bed", sep="\t", index=False, header=False)

    # save split
    for c in CHROMS[:-2]:
        bed[bed["chrom"] == f"chr{c}"].to_csv(f"/mnt/data_3/hallep/nagnag/snps/intersect/canon_chr{c}.bed", sep="\t", index=False, header=False)

    print("\n===== Canonical NAGNAGs =====")
    print(bed)

''' IGV .BED Files '''

# NAME + RGB for 1-OFF NAGNAGs
def bed_1off(site:pd.Series) -> tuple[str, str]:

    ''' Get RGB for 1-off NAGNAGs
    
    Parameters
    ----------
    site : pandas.Series
        1-off splice site information \\
        required column: "index", "off_pos", "motif_seq", "strand", "csite_pos"
    
    Returns
    -------
    name : str
        {index}_off{off pos}_{motif}_{strand}_{splice type}
    _ : str
        {r},{g},{b} values, based on "off_pos"
    '''

    # splice type
    def stype(csites:str):

        # cleavage site positions
        pos = sorted(map(int, csites.split(",")))

        # altenatively-spliced
        if (0 in pos) and (1 in pos):
            return "AS"
        
        # proximally-spliced
        if (0 in pos):
            return "PS"

        # distally-spliced
        if (1 in pos):
            return "DS"
        
        return ""

    name = f"{site["index"]}_off{site["off_pos"]}_{site["motif_seq"].upper()}_{site["strand"]}_{stype(site["csite_pos"])}"

    # NBGNAG: bright orange (xkcd)
    if str(site["off_pos"]) == "1":
        return name, "255,91,0"
    
    # NAHNAG: marigold (xkcd)
    if str(site["off_pos"]) == "2":
        return name, "252,192,6"
    
    # NAGNBG: raspberry (xkcd)
    if str(site["off_pos"]) == "4":
        return name, "176,1,73"

    # NAGNAH: plum (xkcd)
    if str(site["off_pos"]) == "2":
        return name, "8,15,65"

    # default: black
    return name, "0,0,0"

# NAME + RGB for ALL (UPPERCASE) 3' SPLICE SITES
def bed_3ss(site:pd.Series) -> tuple[str, str]:

    ''' Get RGB for 1-off NAGNAGs
    
    Parameters
    ----------
    site : pandas.Series
        1-off splice site information \\
        required column: "index", "ssite_type", "ssite_seq", "strand", "off_pos"
    
    Returns
    -------
    name : str
        {index}_{site type}_{motif}_{strand}
    _ : str
        {r},{g},{b} values, based on "ssite_type"
    '''

    name = f"{site["index"]}_{site["ssite_type"]}_{site["ssite_seq"].upper()}_{site["strand"]}"

    # Canonical 1-NAG: purple
    if str(site["ssite_type"]) == "1C":
        return name, "88,30,49"
    
    # Non-Canonical 1-NAG: orange/brown
    if str(site["ssite_type"]) == "1NC":
        return name, "194,97,0"
    
    # Canonical NAGNAG: red
    if str(site["ssite_type"]) == "2C":
        return name, "165,36,13"

    # Non-Canonical NAGNAG: yellow
    if str(site["ssite_type"]) == "2NC":
        return name, "255,227,59"

    # 3+ NAGs: blue
    if str(site["ssite_type"]) == "3+":
        return name, "50,122,133"

    # default: black
    return name, "0,0,0"

# IGV .BED FILES
def create_igv_bed(src:str, dst:str, start_col:str, end_col:str, func):

    ''' Create .bed file for IVG
    
    Columns: "chrom", "start", "end", "end", "score", "strand", "thickStart", "thickEnd", "rgb"

    Parameters
    ----------
    src : str
        filepath to source .txt file of splice sites \\
        required columns: "chrom", {start_col}, {end_col}, "strand",
            and any columns used by {func}
    dst : str
        filepath to destination .bed file of splice sites
    start_col, end_col : str
        column names for start and end coordinates
    func : function
        function for determining the name and rgb of the splice site
        * parameters: pandas.Series (pandas.DataFrame row)
        * returns: tuple of 2 strings
    '''

    # load splice sites
    txt = pd.read_csv(src, sep="\t", dtype=dtypes)

    names, rgbs = zip(*[func(site[1]) for site in txt.iterrows()])

    # create DataFrame
    bed = pd.DataFrame({
        "chrom" : txt["chrom"].values,
        "start" : txt[start_col].values,
        "end" : txt[end_col].values,
        "name" : names,
        "score" : [0] * len(txt),
        "strand" : txt["strand"].values,
        "thickStart" : txt[start_col].values,
        "thickEnd" : txt[end_col].values,
        "rgb" : rgbs
    })
    
    # save
    bed.to_csv(dst, sep="\t", index=False, header=False)

    print(bed)

''' ===== VCF FILES ===== '''

''' Parse '''

# PARSE .VCF VARIANT
def parse_variant(variant:str) -> dict:

    ''' Parse .VCF variant entry
    
    Parameters
    ----------
    variant : str
        line from .VCF
    
    Returns
    -------
    _ : dict
        variant information
        * keys: "chrom", "pos", "id", "ref", "alt", "qual", "filt", "info"
        * values: information values
    '''

    # commented line
    if variant[0] == "#":
        return None

    return {c : i for c,i in zip(vcf_columns, variant.strip().split("\t"))}

# EXTRACT INFO FIELD VALUE
def get_info_value(info:str, field:str, default:type):
    
    ''' Get value from info column
    
    Parameters
    ----------
    info : str
        entry in VCF "INFO" column
    field : str
        "INFO" field name
    default : type
        default value type for {field}

    Returns
    -------
    _ : any
        info value
    '''

    # no info
    if info == ".":
        return info_defaults[default]
    
    # split fields
    values = {i.split("=")[0] : i.split("=")[1].strip() if ("=" in i) else True for i in info.replace(";_", "%3B_").split(";")}

    # if specified
    if field in values.keys():
        return default(values[field])
    
    # default
    return info_defaults[default]

# PARSE .VCF INFO COLUMN
def parse_info(info:str, info_fields:dict) -> list:
    
    ''' Parse INFO column from .VCF
    
    Parameters
    ----------
    info : str
        value from .VCF
    info_fields : dict
        possible fields in "info"
        * keys: field name
        * values: datatype
    
    Returns
    -------
    _ : list
        values: info field values (see {info_fields} dictionary)
    '''

    # no info
    if info == ".":
        spec = {}

    # split fields
    else:
        spec = {i.split("=")[0] : i.split("=")[1] if ("=" in i) else True for i in info.replace(";_", "%3B_").split(";")}

        # replace special characters    
        for f in ["DNA", "PROT", "PHEN"]:
            if f in spec.keys():
                spec[f] = spec[f].replace("%3A", ":").replace("%3B", ";").replace("%3D", "=").replace("%2C", ",")
                
        # remove extraneous quotation marks
        for f,t in spec.items():
            if t == str:
                spec[f] = spec[f].strip("\"")

    return [spec[f] if f in spec.keys() else info_defaults[t] for f,t in info_fields.items()]

# ADD PARSED 'INFO' COLUMNS
def add_info_columns(df:str|pd.DataFrame, info_fields:dict) -> pd.DataFrame:

    ''' Split "info" column from .vcf into separate columns

    
    Parameters
    ----------
    df : str or pandas.DataFrame
        filepath to .txt from parsed .vcf file OR parsed .vcf file DataFrame \\
        if filepath, saves with new columns \\
        required column: "info"
    info_fields : dict
        possible fields in "info"
        * keys: field name
        * values: datatype
    
    Returns
    -------
    df : pandas.DataFrame
        parsed .vcf file DataFrame with "info" columns
    '''

    save = df if type(df) == str else ""

    # load snps
    if type(df) == str:
        df = pd.read_csv(df, sep="\t", dtype=dtypes)

    # parse "info" column (parallel processing)
    res = zip(*progress_starmap(parse_info, [(l, info_fields) for l in df["info"].values], n_cpu=24))
    
    # add columns
    for col,val in zip(info_fields.keys(), res):
        df[col] = val

    # save
    if save:
        df.to_csv(save, sep="\t", index=False)
        print(df)

    return df

# ADD ALLELE FREQUENCY COLUMNS
def add_frequency_columns(df:str|pd.DataFrame) -> pd.DataFrame:

    ''' Get allele frequencies from dbSNP .vcf

    
    Parameters
    ----------
    df : str or pandas.DataFrame
        filepath to .txt from parsed .vcf file OR parsed .vcf file DataFrame \\
        if filepath, saves with new columns \\
        required columns: "CAF", "TOPMED"
    
    Returns
    -------
    df : pandas.DataFrame
        parsed .vcf DataFrame with allele frequency columns
    '''

    save = df if type(df) == str else ""

    # load snps
    if type(df) == str:
        df = pd.read_csv(df, sep="\t", dtype=dtypes)

    # get frequencies
    df["ref_freq_CAF"] = [float(v.split(",")[0]) if (v != ".") and (v.split(",")[0] != ".") else -1 for v in df["CAF"].values]
    df["alt_freq_CAF"] = [float(v.split(",")[1]) if (v != ".") and (v.split(",")[1] != ".") else -1 for v in df["CAF"].values]
    df["ref_freq_TOPMED"] = [float(v.split(",")[0]) if (v != ".") and (v.split(",")[0] != ".") else -1 for v in df["TOPMED"].values]
    df["alt_freq_TOPMED"] = [float(v.split(",")[1]) if (v != ".") and (v.split(",")[1] != ".") else -1 for v in df["TOPMED"].values]
    
    # save
    if save:
        df.to_csv(save, sep="\t", index=False)
        print(df)

    return df

''' Run '''

# PARSE .VCF FILE
def parse_vcf(src:str, dst:str, info_fields:dict=None, dst_canon:str=None):

    ''' Parse .vcf file into .txt format

    columns: "id", "chrom", "pos", "ref", "alt", "qual", "filt", "info", columns in {info_fields}

    Parameters
    ----------
    src : str
        filepath to source .vcf file
    dst : str
        filepath to destination file of variant information
    info_fields : dict (default = None)
        possible fields in "info" to add to DataFrame \\
        if None, do not parse info columns
        * keys: field name
        * values: datatype
    dst_canon : str (default = None)
        filepath to destination file of variant information for variants on canonical chromosomes \\
        if None, do not check
    '''

    # load VCF
    df = pd.read_csv(src, sep="\t", comment="#", header=None, names=vcf_columns,
                     dtype={0:"str", 1:"int", 2:"str", 3:"str", 4:"str",
                            5:"str", 6:"str", 7:"str"}).set_index(keys="id")

    # parse info columns
    if info_fields:
        df = add_info_columns(df=df, info_fields=info_fields)

    # save
    df.to_csv(dst, sep="\t", index=True)

    # save canonical
    if dst_canon:
        canon += [f"chr{c}" for c in CHROMS]

        cv = df[df["chrom"].apply(lambda x: str(x) in canon)]
        cv.to_csv(dst_canon, sep="\t", index=True)

        print("\n===== Variants on Canonical Chromosomes =====")
        print(cv)
    
    print("\n===== All Variants =====")
    print(df)

# PARSE .VCF FILE BY CHROMOSOME
def parse_vcf_chrom(src:str, dst_chr:str, dst_all:str, n_lines:int=None):

    ''' Parse .vcf file into .txt format
    
    columns: "id", "chrom", "pos", "ref", "alt", "qual", "filt", "info"

    Parameters
    ----------
    src : str
        filepath to source .vcf file
    dst_chr : str
        filepath to destination .txt for variants by chromosome, where chromosome number is indicated by '#'
    dst_all : str
        filepath to destination .txt for all variants
    n_lines : int (default = None)
        total number of lines to parse
    '''
    
    cols = ["id"] + vcf_columns[:2] + vcf_columns[3:]

    # open files
    files = {c : open(dst_chr.replace("#", c), "w") for c in CHROMS[:-1]}
    files["all"] = open(dst_all, "w")

    # data header
    header = "\t".join(cols) + "\n"
    for f in files.values():
        f.write(header)

    # for each .VCF entry:
    for l in tqdm(open(src, "r"), total=n_lines):
        parsed = parse_variant(l)

        # valid variant
        if parsed:
            vl = "\t".join([parsed[c] for c in cols]) + "\n"

            files[parsed["chrom"]].write(vl)
            files["all"].write(vl)

    # close files
    [file.close() for file in files.values()]

''' .bed Files '''

# CREATE VARIANT .BED FILE
def create_vcf_bed(src:str, dst:str, chr_missing:bool=True, use_pd:bool=False):

    ''' Create .bed file for variants
    
    columns: "chrom", "start", "end", "name", "score", "strand"
    * "name" is rs (from "id" column)
    * "score" is the alternate base(s)
    * "strand" is "."

    Parameters
    ----------
    src : str
        path to source .txt file
    dst : str
        path to destination .bed file
    chr_missing : bool (default = True)
        whether the prefix 'chr' is missing from "chrom" column
    use_pd : bool (default = False)
        whether to use pandas.DataFrame to create .bed file
    '''
    
    # create with pandas.DataFrame
    if use_pd:

        # load .txt file
        txt = pd.read_csv(src, sep="\t", dtype=dtypes)

        # format columns
        if chr_missing:
            txt["chrom"] = "chr" + txt["chrom"]
        txt["start"] = txt["pos"]-1
        txt["strand"] = ["."] * len(txt)

        # save as .bed
        txt.to_csv(dst, sep="\t", index=False, header=False,
                   columns=["chrom", "start", "pos", "id", "alt", "strand"])

        print(txt)

    # parse line by line
    else:
        with open(dst, "w") as bed:

            # for each line:
            for i,line in tqdm(enumerate(open(src, "r"))):

                # if not header:
                if i > 0:
                    l = line.strip().split("\t")
                    rs, chrom, pos, alt = l[0], l[1], l[2], l[4]
                    bed.write(f"{"chr" if chr_missing else ""}{chrom}\t{int(pos)-1}\t{pos}\t{rs}\t{alt}\t{"."}\n")

# CREATE VARIANT .BED FILES FOR ALL CHROMOSOMES
def create_vcf_chrom_bed(src:str, dst:str, chr_missing:bool=True):

    ''' Create .BED file for variants on each chromosome
    
    call create_vcf_bed()

    Parameters
    ----------
    src : str
        filepath to .txt of variants by chromosome, where chromosome number is indicated by '#'
    dst : str
        filepath to destination .bed for variants by chromosome, where chromosome number is indicated by '#'
    chr_missing : bool (default = True)
        whether the prefix 'chr' is missing from "chrom" column
    '''

    # create chromosome .bed files
    for c in CHROMS[:-1]:
        create_vcf_bed(src.replace("#", c), dst.replace("#", c), chr_missing=chr_missing, use_pd=False)

''' ===== SQL FILES ===== '''

# EXPORT MYSQL TABLES
def export_database(database:str, files:str|list, tables:str|list=None, host:str="Dexter", user:str="root"):

    ''' export MySQL table to tab-separated .txt files
    
    Password for MySQL session is prompted in the terminal

    Parameters
    ----------
    database : str
        database in which table is stored
    files : str (default = None)
        filepath(s) to destination tab-separated .txt file(s) \\
        if multiple tables and 1 filename, replace # with table name
    tables : str or list (default = None)
        name of table(s) to export \\
        if None, export all tables
    host : str (default = "Dexter")
        host name
    user : str (default = "root")
        user name
    '''

    # get password
    password = maskpass.askpass(mask="")

    # connect to database
    engine = create_engine(f"mysql+mysqlconnector://{user}:{password}@{host}/{database}")

    # export all tables
    if tables == None:
        metadata = MetaData()
        metadata.reflect(engine)
        tables = list(metadata.tables.keys())

    # define {file} names
    if isinstance(tables, list) and isinstance(files, str):
        files = [files.replace("#", t) for t in tables]

    # format {tables} and {files} as lists
    if isinstance(tables, str):
        tables = [tables]
    if isinstance(files, str):
        files = [files]
        
    # for each table:
    for t,f in zip(tables, files):

        # export
        txt = pd.read_sql(sql=f"SELECT * FROM {t}", con=engine)
        txt.to_csv(f, sep="\t", index=False)

        print(f"===== {t} =====")
        print(txt)

# PROCESS HGMD TABLE (.TXT, .BED)
def process_sql(src:str, dst_txt:str, dst_bed:str, id_col:str="acc_num", drop_cols:list[str]=[]):

    ''' Process HGMD table

    source files:
    * variants: {src}
    * coordinates: "/mnt/data_3/hallep/reference/snp_data/HGMD/HGMD_pro/hg38_coords.txt"

    destination files:
    * {dst_txt}
        * only variants (rows) with coordinates (in hg38_coords.txt) are kept
        * index is {id_col}
        * columns:
            * {src} columns: all except {drop_cols}
            * new columns: "chrom", "start", "end", "strand", "ref", "alt"
    * {dst_bed}
        * columns: "chrom", "start", "end", "name", "score", "strand"
            * "name" is id (from {id_col}/"acc_num")
            * "score" is alternative allele

    Parameters
    ----------
    src : str
        filepath to source .txt file containing HGMD table
    dst_txt : str
        filepath to destination .txt file to which to save processed table
    dst_bed : str
        filepath to destination .bed file for variants in table
    id_col : str (default = "acc_num")
        column in {src} containing variant id that matches "acc_num" in hg38_coords.txt
    drop_cols : list of str (default = [])
        columns in {src} to remove from {dst_txt}
    '''
    
    # load coordinates
    hg38 = pd.read_csv("/mnt/data_3/hallep/reference/snp_data/HGMD/HGMD_pro/hg38_coords.txt", sep="\t", index_col="acc_num")
    hg38["chromosome"] = "chr" + hg38["chromosome"]

    # load table, drop columns, remove variants without coordinates
    table = pd.read_csv(src, sep="\t")
    table.drop(columns=drop_cols, inplace=True)
    table = table[table[id_col].apply(lambda x: x in hg38.index)]

    # get coordinates
    table["chrom"] = [hg38.loc[i]["chromosome"] for i in table[id_col]]
    table["start"] = [hg38.loc[i]["coordSTART"] for i in table[id_col]]
    table["end"] = [hg38.loc[i]["coordEND"] for i in table[id_col]]
    table["strand"] = [hg38.loc[i]["strand"] for i in table[id_col]]
    table["ref"] = table["base"].apply(lambda x: x.split("-" if "-" in x else ">")[0])
    table["alt"] = table["base"].apply(lambda x: x.split("-" if "-" in x else ">")[1])

    print(table[table["strand"] == "+"])
    print(table[table["strand"] == "-"])

    # save .bed file
    bed = table.copy()
    bed["start"] = bed["start"] - 1

    bed.to_csv(dst_bed, sep="\t", columns=["chrom", "start", "end", id_col, "alt", "strand"], index=False, header=False)

    print(f"\n===== {dst_bed.split("/")[-1]} =====")
    print(bed[["chrom", "start", "end", id_col, "alt", "strand"]])

    # save valid variants
    table.set_index(keys=id_col, drop=True, inplace=True)
    table.to_csv(dst_txt, sep="\t", index=True)

    print(f"\n===== {dst_txt.split("/")[-1]} =====")
    print(table)

# ADD INFORMATION COLUMNS
def add_allmut_columns(df:str|pd.DataFrame):

    ''' Split "info" column from .vcf into separate columns

    Note: either {file} or {df} must be specified
    
    Parameters
    ----------
    df : str or pandas.DataFrame
        filepath to .txt from parsed SQL table OR parsed SQL table DataFrame
    
    Returns
    -------
    df : pandas.DataFrame
        parsed .vcf file DataFrame with "info" columns
    '''

    columns = ["gnomad_AF", "rankscore"]

    # load allmut table
    allmut = pd.read_csv("/mnt/data_3/hallep/reference/snp_data/HGMD/HGMD_pro/allmut.txt", sep="\t",
                         index_col="acc_num", usecols=["acc_num"] + columns)

    save = df if type(df) == str else ""

    # load snps
    if type(df) == str:
        df = pd.read_csv(df, sep="\t", dtype=dtypes)

    # add columns
    for col in columns:
        df[col] = [allmut.loc[a][col] for a in df["acc_num"].values]

    # save
    if save:
        df.to_csv(save, sep="\t", index=False)
        print(df)

    return df

''' ===== RUN ===== '''

''' 3' splice sites

create_3ss_bed()
create_nagnag_1off_bed()

create_igv_bed(src="/mnt/data_3/hallep/nagnag/snps/1off_splice_sites.txt",
               dst="/mnt/data_3/hallep/reference/tx_data/1off_nagnags.bed",
               start_col="motif_start", end_col="motif_end", func=bed_1off)

create_igv_bed(src="/mnt/data_3/hallep/nagnag/3ss_uppercase.txt",
               dst="/mnt/data_3/hallep/reference/tx_data/3ss_uppercase.bed",
               start_col="ssite_start", end_col="ssite_end", func=bed_3ss)
# '''

''' dbSNP: VCF

parse_vcf_chrom(src="/mnt/data_3/hallep/reference/snp_data/dbSNP/00-All.vcf",
                dst_chr="/mnt/data_3/hallep/reference/snp_data/dbSNP/chrom_split/dbSNP_chr#.txt",
                dst_all="/mnt/data_3/hallep/reference/snp_data/dbSNP/all_dbSNP_snps.txt", n_lines=660146231)

create_vcf_chrom_bed(src="/mnt/data_3/hallep/reference/snp_data/dbSNP/chrom_split/dbSNP_chr#.txt",
                     dst="/mnt/data_3/hallep/reference/snp_data/bed/dbsnp_chr#.bed", chr_missing=True)
# '''

''' HGMD (VCF)

parse_vcf(src="/mnt/data_3/hallep/reference/snp_data/HGMD/HGMD_Bundle_2023/hgmd_pro_2023.1_hg38.vcf",
          dst="/mnt/data_3/hallep/reference/snp_data/HGMD/all_HGMD_snps.txt", info_fields=info_fields["HGMD"])

parse_vcf(src="/mnt/data_3/hallep/reference/snp_data/HGMD/HGMD_Bundle_2023/hgmd_pro_2023.1_hg19.vcf",
          dst="/mnt/data_3/hallep/reference/snp_data/HGMD/all_HGMD_snps_hg19.txt", info_fields=info_fields["HGMD"])

create_vcf_bed(src="/mnt/data_3/hallep/reference/snp_data/HGMD/all_HGMD_snps.txt",
               dst="/mnt/data_3/hallep/reference/snp_data/bed/hgmd_vcf.bed", chr_missing=True, use_pd=True)
# '''

''' HGMD Splice (SQL)

# export_database(database="hgmd_pro", tables=None,
#                 files=f"/mnt/data_3/hallep/reference/snp_data/HGMD/HGMD_pro/#.txt",
#                 host="localhost", user="root")

# export_database(database="hgmd_phenbase", tables=None,
#                 files=f"/mnt/data_3/hallep/reference/snp_data/HGMD/HGMD_phenbase/#.txt",
#                 host="localhost", user="root")

# export_database(database="hgmd_snp", tables=None,
#                 files=f"/mnt/data_3/hallep/reference/snp_data/HGMD/HGMD_snp/#.txt",
#                 host="localhost", user="root")

# export_database(database="hgmd_views", tables=None,
#                 files=f"/mnt/data_3/hallep/reference/snp_data/HGMD/HGMD_views/#.txt",
#                 host="localhost", user="root")

process_sql(src="/mnt/data_3/hallep/reference/snp_data/HGMD/HGMD_pro/splice.txt",
            dst_txt="/mnt/data_3/hallep/reference/snp_data/HGMD/hgmd_splice_variants.txt",
            dst_bed="/mnt/data_3/hallep/reference/snp_data/bed/hgmd_splice.bed",
            drop_cols=["author", "journal", "fullname", "vol", "page", "year", "pmid"])
# '''

''' ClinVar (VCF)

parse_vcf(src="/mnt/data_3/hallep/reference/snp_data/ClinVar/clinvar.vcf",
          dst="/mnt/data_3/hallep/reference/snp_data/ClinVar/all_clinvar_snps.txt", info_fields=None,
          dst_canon="/mnt/data_3/hallep/reference/snp_data/ClinVar/canon_clinvar_snps.txt")

create_vcf_bed(src="/mnt/data_3/hallep/reference/snp_data/ClinVar/canon_clinvar_snps.txt",
               dst="/mnt/data_3/hallep/reference/snp_data/bed/clinvar.bed", chr_missing=True, use_pd=True)
# '''
