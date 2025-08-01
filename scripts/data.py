''' Download NCBI RefSeq transcript data and hg38 reference genome sequences '''

import io
import requests
import gzip
import pickle
from Bio import SeqIO

# Download NCBI RefSeq data
def download_ncbi_refseq():

    ''' Download NCBI RefSeq data
    
    src: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz
    
    dst: source/ncbiRefSeq.txt    
    '''

    # download file
    gz_file = requests.get("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz")

    # extract content
    with gzip.GzipFile(fileobj=io.BytesIO(gz_file.content)) as gz_content:
        content = gz_content.read()

    # write to file
    with open("source/ncbiRefSeq.txt", "w") as file:
        file.write(content.decode('utf-8'))

# Download hg38 chromosome sequences
def download_hg38():

    ''' Download hg38 sequence fasta file
    
    src: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    
    dst: source/hg38.fa
    '''

    # download file
    gz_file = requests.get("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz")

    # extract content
    with gzip.GzipFile(fileobj=io.BytesIO(gz_file.content)) as gz_content:
        content = gz_content.read()

    # write to file
    with open("source/hg38.fa", "w") as file:
        file.write(content.decode('utf-8'))

# Create sequence dictionary
def create_hg38_seq_dict():
    
    ''' Create a sequence dictionary for reference human genome hg38
    
    src: source/hg38.fa
    
    dst: source/hg38 (pickle)

        keys: str
            chromosome name: chr#
        values: Bio.Seq.Seq
            0-indexed chromosome sequence  
    '''

    hg38_dict = {}

    # create dictionary
    for record in SeqIO.parse("source/hg38.fa", "fasta"):
        hg38_dict[record.id] = record.seq
    
    # pickle
    with open("source/hg38", "wb") as file:
        pickle.dump(hg38_dict, file)
