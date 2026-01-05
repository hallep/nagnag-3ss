#!/bin/bash

print_help() {
    echo -e "Usage: run.sh [-h] [-t INT] [-s FILE] [-c FILE] [-x] [-H]\n"
    echo -e "Run NAGNAG 3' splice site analyses \n"

    echo -e "Options:"
    echo -e "   -h\t\tShow this message and exit"
    echo -e "   -t INT\tNumber of parallel threads available for multiprocessing [24]"
    echo -e "   -s FILE\tPath to tab-delimited HGMD splice file [$PWD/src/HGMD_pro/splice.txt]"
    echo -e "   -c FILE\tPath to tab-deliimted HGMD hg38_coords file [$PWD/src/HGMD_pro/hg38_coords.txt]"
    echo -e "   -H\t\tDo not process or analyze HGMD Splice variants"

    echo -e "\nTo generate HGMD files, run 'export_hgmd.py -u <user> -h <host>'"

    exit 0
}

if [[ $1 == "-h" || $1 == "--help" ]]; then
    print_help
fi

# HGMD variant files
while getopts "h:t:s:c:H" opt; do
    case $opt in
        h) print_help ;;
        t) THREADS="--num-threads $OPTARG" ;;
        s) SPLICE="--splice $OPTARG" ;;
        c) COORDS="--coords $OPTARG" ;;
        H) HARG="-H" ;;
    esac
done

(
echo "[$(date +"%Y-%m-%d %H:%M:%S")] CMD: run.sh $@" >&2

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
python nagnag/09-find-variants.py $THREADS $HARG

mkdir -p variants/stats
python nagnag/10-count-variants.py $THREADS $HARG

mkdir -p variants/proteome
python nagnag/11-analyze-variants.py $HARG

# Done
echo "[$(date +"%Y-%m-%d %H:%M:%S")] ---------- DONE ----------" >&2

) 2>&1 | tee nagnag.log
