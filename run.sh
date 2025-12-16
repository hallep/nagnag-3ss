#!/bin/bash

echo "[$(date +"%Y-%m-%d %H:%M:%S")] CMD: run.sh $@" >&2

# HGMD variant files
while getopts "s:c:H" opt; do
    case $opt in
        s) SPLICE="--splice $OPTARG" ;;
        c) COORDS="--coords $OPTARG" ;;
        H) HARG="-H"
    esac
done

# Download data
mkdir -p src
python nagnag/00-data.py $HARG $SPLICE $COORDS

# Identify splice sites
mkdir -p sites
python nagnag/01-sites.py

# Calculate frequency statistics
mkdir -p stats
python nagnag/02-stats.py

# Create motif heatmaps + logos
mkdir -p figures/heatmaps
python nagnag/03-motifs.py

mkdir -p figures/logos
python nagnag/04-logos.py

# Analyze proteomic effects
mkdir -p proteome
python nagnag/05-proteome.py

python nagnag/06-poswise.py
python nagnag/07-splice-selection.py

# Variants
mkdir -p variants/src variants/bed
python nagnag/08-init-variants.py

mkdir -p variants/affecting variants/found
python nagnag/09-find-variants.py $HARG

mkdir -p variants/stats
python nagnag/10-count-variants.py $HARG

mkdir -p variants/proteome
python nagnag/11-analyze-variants.py $HARG

# Done
echo "[$(date +"%Y-%m-%d %H:%M:%S")] ---------- DONE ----------" >&2
