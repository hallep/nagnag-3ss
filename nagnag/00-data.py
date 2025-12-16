''' Download data

**Created Files:**
* `src/ncbiRefSeq.txt`
* `src/hg38.fa`
* `src/hg38`
* `src/dbSNP.vcf`
* `src/ClinVar.vcf`
'''

from lib import argparse, subprocess, pickle, Path, SeqIO, tqdm, pd
from utils import ROOT, log_script, log_fn
from utils.seq import CHROMS
from utils.vnt import txt_cols, vcf_cols

# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-H", "--ignore-HGMD", action="store_true", help="Do not process/analyze HGMD Splice data")
parser.add_argument("-s", "--splice", metavar="FILE", type=Path, default=Path(f"{ROOT}/src/HGMD_pro/splice.txt"),
                    help=f"path to HGMD splice file [{ROOT}/src/HGMD_pro/splice.txt]")
parser.add_argument("-c", "--coords", metavar="FILE", type=Path, default=Path(f"{ROOT}/src/HGMD_pro/hg38_coords.txt"),
                    help=f"path to HGMD hg38_coords file [{ROOT}/src/HGMD_pro/hg38_coords.txt]")
args = parser.parse_args()

# get data
def download(name:str, link:str):
    subprocess.run(["wget", "-O", f"{ROOT}/src/{name}", link])
    subprocess.run(["gunzip", "-f", f"{ROOT}/src/{name}"])

# ----- Process Variants ----- #

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

# ===== RUN ===== #
log_script("00-data.py")

# NCBI RefSeq
log_fn("Downloading NCBI RefSeq transcripts")
download("ncbiRefSeq.txt.gz", "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz")

# hg38
log_fn("Downloading GRCh38/hg38 reference genome")
download("hg38.fa.gz", "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz")

log_fn("Parsing GRCh38/hg38 .fasta file")
hg38_dict = {record.id : record.seq for record in SeqIO.parse(f"{ROOT}/src/hg38.fa", "fasta")}
pickle.dump(hg38_dict, (open(f"{ROOT}/src/hg38", "wb")))

# dbSNP
log_fn("Downloading dbSNP variants")
download("dbSNP.vcf.gz", "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz")

log_fn("Processing dbSNP .vcf file")
vcf2txt2bed(vcf=f"{ROOT}/src/dbSNP.vcf", txt=f"{ROOT}/variants/src/dbSNP_$_#.txt",
            bed=f"{ROOT}/variants/bed/dbSNP_$_#.bed", by_chrom=True, split_freq=True)

# ClinVar
log_fn("Downloading ClinVar variants")
download("ClinVar.vcf.gz", "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz")

log_fn("Processing ClinVar .vcf file")
vcf2txt2bed(vcf=f"{ROOT}/src/ClinVar.vcf", txt=f"{ROOT}/variants/src/ClinVar.txt",
            bed=f"{ROOT}/variants/bed/ClinVar.bed", by_chrom=False, split_freq=False)

# HGMD Splice
if not args.ignore_hgmd:
    log_fn("Processing HGMD Splice SQL data")
    sql2txt2bed(sql=args.splice, txt=f"{ROOT}/variants/src/HGMD_splice.txt",
                bed=f"{ROOT}/variants/bed/HGMD_splice.bed", coords=args.coords)
