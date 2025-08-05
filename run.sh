#!/bin/bash

# ===== GET SOURCE DATA ===== #
mkdir -p src

# NCBI RefSeq Data
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz
gunzip ncbiRefSeq.txt.gz
mv ncbiRefSeq.txt src

# Chromosome Sequences
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
mv hg38.fa src
python scripts/data.py

# dbSNP Variants
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz
gunzip 00-All.vcf.gz
mv 00-All.vcf src/dbSNP.vcf

# ClinVar Variants
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
gunzip clinvar.vcf.gz
mv clinvar.vcf src

# ===== IDENTIFY SPLICE SITES ===== #
mkdir -p sites
python scripts/sites.py

# ===== ANALYZE PROTEOME EFFECTS ===== #
mkdir -p stats
python scripts/stats.py

# ===== ANALYZE PROTEOME EFFECTS ===== #
python scripts/proteome.py

# ===== CONTROL FOR SPLICE SELECTION ===== #
mkdir -p proteome
python scripts/poswise.py
python scripts/splice_selection.py
