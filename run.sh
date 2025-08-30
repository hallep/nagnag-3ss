#!/bin/bash

# ===== GET SOURCE DATA ===== #
mkdir -p src

# NCBI RefSeq Data
wget -O src/ncbiRefSeq.txt.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz
gunzip src/ncbiRefSeq.txt.gz
echo "NCBI RefSeq transcripts downloaded"

# Chromosome Sequences
wget -O src/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip src/hg38.fa.gz
echo "GRCh38/hg38 chromosome sequences downloaded"
python scripts/data.py

# dbSNP Variants
wget -O src/dbSNP.vcf.gz https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz
gunzip src/dbSNP.vcf.gz
echo "dbSNP variants downloaded"

# ClinVar Variants
wget -O src/ClinVar.vcf.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
gunzip src/ClinVar.vcf.gz
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
mkdir -p variants/src && mkdir -p variants/bed
python scripts/init_variants.py

mkdir -p variants/affecting && mkdir -p variants/found
python scripts/find_variants.py

mkdir -p variants/stats
python scripts/count_variants.py

mkdir -p variants/proteome
python scripts/analyze_variants.py
