''' Download sequence, transcript, and variant data '''

from lib import ROOT, requests, io, gzip, pickle, SeqIO

# NCBI RefSeq data
def ncbi_refseq():

    ''' Download NCBI RefSeq transcript annotation data
    
    Source: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz
    
    Destination: `src/ncbiRefSeq.txt`
    '''

    # download file
    gz_file = requests.get("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz")
    
    # unzip and extract content
    with gzip.GzipFile(fileobj=io.BytesIO(gz_file.content)) as gz_content:
        content = gz_content.read()

    # save
    with open(f"{ROOT}/src/ncbiRefSeq.txt", "w") as txt_file:
        txt_file.write(content.decode())

# GRCh38/hg38 chromosome sequences
def hg38_sequence():
    
    ''' Download GRCh38/hg38 chromosome sequences and save as dictionary
    
    Source: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    
    Destination:
    * fasta: `src/hg38.fa`
    * dictionary: `src/hg38`
    '''

    # download file
    gz_file = requests.get("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz")

    # unzip and extract content
    with gzip.GzipFile(fileobj=io.BytesIO(gz_file.content)) as gz_content:
        content = gz_content.read()

    # save .fasta
    with open(f"{ROOT}/src/hg38.fa", "w") as fa_file:
        fa_file.write(content.decode())

    # convert to dictionary
    hg38_dict = {record.id : record.seq for record in SeqIO.parse(f"{ROOT}/src/hg38.fa", "fasta")}

    # save dictionary
    pickle.dump(hg38_dict, (open(f"{ROOT}/src/hg38", "wb")))

# dbSNP variants
def dbsnp_variants():

    ''' Download dbSNP variants 
    
    Source: https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz

    Destination: `src/dbSNP_variants.vcf`
    '''

    # download file
    gz_file = requests.get("https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz")
    
    # unzip and extract content
    with gzip.GzipFile(fileobj=io.BytesIO(gz_file.content)) as gz_content:
        content = gz_content.read()

    # save
    with open(f"{ROOT}/src/dbSNP_variants.vcf", "w") as vcf_file:
        vcf_file.write(content.decode())

# ClinVar variants
def clinvar_variants():

    ''' Download ClinVar variants 
    
    Source: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

    Destination: `src/ClinVar_variants.vcf`
    '''

    # download file
    gz_file = requests.get("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz")
    
    # unzip and extract content
    with gzip.GzipFile(fileobj=io.BytesIO(gz_file.content)) as gz_content:
        content = gz_content.read()

    # save
    with open(f"{ROOT}/src/ClinVar_variants.vcf", "w") as vcf_file:
        vcf_file.write(content.decode())

ncbi_refseq()
hg38_sequence()
dbsnp_variants()
clinvar_variants()
