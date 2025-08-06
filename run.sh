#!/bin/bash

# ===== GET SOURCE DATA ===== #
mkdir -p src

# NCBI RefSeq Data
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz
gunzip ncbiRefSeq.txt.gz
mv ncbiRefSeq.txt src
echo "NCBI RefSeq transcripts downloaded"

# Chromosome Sequences
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
mv hg38.fa src
echo "GRCh38/hg38 chromosome sequences downloaded"
python scripts/data.py

# dbSNP Variants
wget https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz
gunzip 00-All.vcf.gz
mv 00-All.vcf src/dbSNP.vcf
echo "dbSNP variants downloaded"

# ClinVar Variants
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
gunzip clinvar.vcf.gz
mv clinvar.vcf src
echo "ClinVar variants downloaded"

# ===== IDENTIFY SPLICE SITES ===== #
mkdir -p sites
python scripts/sites.py

# ===== RUN FREQUENCY STATISTICS ===== #
mkdir -p stats
python scripts/stats.py

# ===== CREATE MOTIF HEATMAPS + LOGOS ===== #
mkdir -p figures/heatmaps
python scripts/motifs.py

mkdir -p figures/logos
python scripts/logos.py

# ===== ANALYZE PROTEOME EFFECTS ===== #
python scripts/proteome.py

# ===== CONTROL FOR SPLICE SELECTION ===== #
mkdir -p proteome
python scripts/poswise.py
python scripts/splice_selection.py

# ===== Variants ===== #
mkdir -p variants/bed && mkdir -p variants/src
python scripts/init_variants.py

mkdir -p variants/found
python scripts/find_variants.py
