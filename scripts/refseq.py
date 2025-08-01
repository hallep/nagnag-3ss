''' Parse NCBI RefSeq transcript annotation data;
extract cleavage sites, splice sites, and splice scenarios '''

import pandas as pd
from scripts.sequence import get_3cs_tri, get_ss_seq

# Get cleavage sites
def get_cleavage_sites():

    ''' Identify all unique 3' cleavage sites
    
    src: source/ncbiRefSeq.tx    
    
    dst: data/3cs.txt

        "index": int
            0-indexed row number in DataFrame
        "chrom": str
            chromosome on which 3' cleavage site is found
            format: chr#
        "down_start": int
            0-indexed chromosome position of start of downstream exon
             * plus strand: first base of downstream exon 
             * minus strand: last base of upstream intron
        "strand": str ("+", "-")
            strand of transcript in which cleavage site is found
        "num_scens": int
            number of splice scenarios associated with each cleavage site
        "rtype": str {"CDS", "5'UTR", "3'UTR", "ncRNA"}
            comma-separated list, 1 for each splice scenario, of the RNA type of cleavage site and transcript
             * "CDS: found in a coding (NM, XM) transcript, within the coding sequence
             * "5'UTR": found in a coding (NM, XM) transcript, upstream of the coding sequence
             * "3'UTR": found in a coding (NM, XM) transcript, downstream of the coding sequence
             * "ncRNA":  found in an non-coding (NR, XR) transcript
        "up_end": str
            comma-separate list, 1 for each splice scenario, of the 0-indexed chromosome position of end of upstream exon
             * plus strand: first base of downstream intron 
             * minus strand: last base of upstream exon
        "frame": str
            comma-separated list, 1 for each splice scenario, of the frame/phase of downstream exon
             * -1: downstream exon is in a UTR or in a non-coding RNA transcript
             * 0: downstream exon is the first base of its codon
             * 1: downstream exon is the second base of its codon
             * 2: downstream exon is the third base of its codon
        "accession": str
            comma-separated list, 1 for each splice scenario, of forward-slash-separated lists of accession numbers
        "n_isoforms": str
            comma-separated list, 1 for each splice scenario, of the number of transcripts associated with each splice scenario    
    '''

    # load transcript data table
    columns = ["bin", "name", "chrom", "strand", "txStart", "txEnd",
               "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
               "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"]
    ncbi_table = pd.read_csv("source/ncbiRefSeq.txt", sep="\t", header=None, names=columns)

    canon_chroms = [f"chr{c}" for c in list(range(1,23)) + ["X", "Y", "M"]]

    # get transcripts on canonical chromosomes
    ncbi_table = ncbi_table.loc[ncbi_table["chrom"].isin(canon_chroms)]
    ncbi_table.reset_index(inplace=True)

    cs3 = {}

    # for each transcript:
    for _,tx in ncbi_table.iterrows():

        # split exon definitions
        exonStarts = [int(s) for s in tx["exonStarts"].strip(",").split(",")]
        exonEnds = [int(e) for e in tx["exonEnds"].strip(",").split(",")]
        exonFrames = tx["exonFrames"].strip(",").split(",")

        # plus (+) strand:
        if tx["strand"] == "+":

            # for each 3' splice site (start of all but first exon):
            for j in range(1, tx["exonCount"]):

                # csite (cleavage site) = (chrom, down_start, strand)
                csite = (tx["chrom"], exonStarts[j], tx["strand"])

                # rtype = RNA Type {"CDS", "5'UTR", "3'UTR", "ncRNA"}
                # non-coding RNA
                if tx["name"][:2] == "XR" or tx["name"][:2] == "NR":
                    rtype = "ncRNA"
                # coding RNA
                else:
                    # coding sequence
                    if exonStarts[j] >= tx["cdsStart"] and exonStarts[j] <= tx["cdsEnd"]:
                        rtype = "CDS"
                    # 5' UTR
                    elif exonStarts[j] < tx["cdsStart"]:
                        rtype = "5'UTR"
                        exonFrames[j] = "-1"
                    # 3' UTR
                    elif exonStarts[j] > tx["cdsEnd"]:
                        rtype = "3'UTR"
                        exonFrames[j] = "-1"
                    
                # scen (splice scenario) = (rtype, down_frame, up_end)
                scen = (rtype, exonFrames[j], exonEnds[j-1])

                # if new cleavage site:
                if csite not in cs3:
                    cs3[csite] = {}
                
                # if new splice scenario:
                if scen not in cs3[csite]:
                    cs3[csite][scen] = list()
                
                # add accession number
                cs3[csite][scen].append(tx["name"])
        
        # minus (-) strand:
        elif tx["strand"] == "-":

            # for each 3' splice site (end of all but last exon)
            for j in range(tx["exonCount"]-1):

                # csite (cleavage site) = (chrom, down_start, strand)
                csite = (tx["chrom"], exonEnds[j], tx["strand"])

                # rtype = RNA Type {"CDS", "5'UTR", "3'UTR", "ncRNA"}
                # non-coding RNA
                if tx["name"][:2] == "XR" or tx["name"][:2] == "NR":
                    rtype = "ncRNA"
                # coding RNA
                else:
                    # coding sequence
                    if exonEnds[j] >= tx["cdsStart"] and exonEnds[j] <= tx["cdsEnd"]:
                        rtype = "CDS"
                    # 3' UTR
                    elif exonEnds[j] < tx["cdsStart"]:
                        rtype = "3'UTR"
                        exonFrames[j] = "-1"
                    # 5' UTR
                    elif exonEnds[j] > tx["cdsEnd"]:
                        rtype = "5'UTR"
                        exonFrames[j] = "-1"

                # scen (splice scenario) = (rtype, down_frame, up_end)
                scen = (rtype, exonFrames[j], exonStarts[j+1])

                # if new cleavage site:
                if csite not in cs3:
                    cs3[csite] = {}
                
                # if new splice scenario:
                if scen not in cs3[csite]:
                    cs3[csite][scen] = list()
                
                # add accession number
                cs3[csite][scen].append(tx["name"])
    
    # parse cleavage sites
    chrom = []
    down_start = []
    strand = []
    num_scens = []
    rtype = []
    up_end = []
    frame = []
    accession = []
    n_isoforms = []

    # for each csite
    for csite in cs3.keys():
        chrom.append(csite[0])
        down_start.append(csite[1])
        strand.append(csite[2])

        num_scens.append(len(cs3[csite]))

        rtype.append(",".join([scen[0] for scen in cs3[csite].keys()]))
        frame.append(",".join([scen[1] for scen in cs3[csite].keys()]))
        up_end.append(",".join([str(scen[2]) for scen in cs3[csite].keys()]))

        accession.append(",".join(["/".join(cs3[csite][scen]) for scen in cs3[csite].keys()]))
        n_isoforms.append(",".join([str(len(cs3[csite][scen])) for scen in cs3[csite].keys()]))

    # create DataFrame    
    csites = pd.DataFrame({
        "chrom" : chrom,
        "down_start" : down_start,
        "strand" : strand,
        "num_scens" : num_scens,
        "rtype" : rtype,
        "up_end" : up_end,
        "frame" : frame,
        "accession" : accession,
        "n_isoforms" : n_isoforms
    })

    # add "index" column
    csites = csites.rename_axis("index").reset_index()

    # save DataFrame
    csites.to_csv("data/3cs.txt", sep="\t", index=False)
    
# Add cleavage site sequences
def add_csite_sequences():

    ''' Add sequences to 3' cleavage site table
    
    file: data/3cs.txt
    
        "csite_seq": str
            last 3 bases of intron
        "scen_eflank_5ss": str
            comma-separated list, 1 for each splice scenario, of the last 100 bases of the upstream exon
        "scen_iflank_5ss": str
            comma-separated list, 1 for each splice scenario, of the first 100 bases of the intron following the upstream exon
        "csite_iflank_3ss": str
            last 100 bases of the intron (including csite_seqs)
        "csite_eflank_3ss": str
            first 100 bases of the downstream exon
    '''

    # load cleavage sites
    csites = pd.read_csv("data/3cs.txt", sep="\t")

    csite_seq = [""] * len(csites)
    eflank_5ss = [""] * len(csites)
    iflank_5ss = [""] * len(csites)
    iflank_3ss = [""] * len(csites)
    eflank_3ss = [""] * len(csites)

    # for each cleavage site:
    for i,cs in csites.iterrows():

        # 3' cleavage site sequences
        csite_seq[i] = get_3cs_tri(cs["chrom"], cs["down_start"], cs["strand"])
        iflank_3ss[i] = get_ss_seq(cs["chrom"], cs["down_start"], cs["strand"], 100, 0)
        eflank_3ss[i] = get_ss_seq(cs["chrom"], cs["down_start"], cs["strand"], 0, 100)

        # 5' cleavage site sequences
        fe = [""] * cs["num_scens"]
        fi = [""] * cs["num_scens"]

        # for each upstream exon end:
        for j,pos in enumerate(cs["up_end"].split(",")):
            fe[j] = get_ss_seq(cs["chrom"], int(pos), cs["strand"], 100, 0)
            fi[j] = get_ss_seq(cs["chrom"], int(pos), cs["strand"], 3, 100)

        eflank_5ss[i] = ",".join(fe)
        iflank_5ss[i] = ",".join(fi)

    # add columns to DataFrame
    csites["csite_seq"] = csite_seq
    csites["scen_eflank_5ss"] = eflank_5ss
    csites["scen_iflank_5ss"] = iflank_5ss
    csites["csite_iflank_3ss"] = iflank_3ss
    csites["csite_eflank_3ss"] = eflank_3ss

    # save DataFrame
    csites.to_csv("data/3cs.txt", sep="\t", index=False)

# Get splice site type
def get_sstype(ssite_seq:str) -> str:

    ''' Get splice site type
    
    Parameters
    ----------
    ssite_seq : str
        sequence of splice site
    
    Returns
    -------
    _ : str {"1C", "1NC", "2C", "2NC", "3+"}
        splice site type
         * "1C": 1-tri, canonical (NAG)
         * "1NC": 1-tri, non-canonical (NBG, NAH, NBH)
         * "2C": 2-tri, canonical (NAGNAG)
         * "2NC": 2-tri, non-canonical
         * "3+": 3+ trinucleotides
    '''

    # 1-tri
    if len(ssite_seq) == 3:
        # canonical
        if (ssite_seq[1].upper() == "A") and (ssite_seq[2].upper() == "G"):
            return "1C"

        # non-canonical
        return "1NC"

    # 2-tri
    elif len(ssite_seq) == 6:

        # canonical
        if ((ssite_seq[1].upper() == "A") and (ssite_seq[2].upper() == "G") and
            (ssite_seq[4].upper() == "A") and (ssite_seq[5].upper() == "G")):
            return "2C"
        
        # non-canonical
        return "2NC"

    # 3+ tri
    return "3+"

# Group cleavage sites
def group_csites():

    ''' Group 3' cleavage sites into 3' splice sites
    
    src: data/3cs.txt
    
    dst: data/3ss.txt
        
        "index": int
            0-indexed row number in DataFrame
        "chrom": str
            chromosome on which 3' splice site is found
            format: chr#
        "strand": str ("+", "-")
            strand of transcript in which splice site is found
        "num_tri": int
            number of cleavage sites + adjacent NAGs in the splice site
        "num_csites": int
            number of cleavage sites in the splice site
        "csite_starts": str
            comma-separated list of cleavage site down_start positions
        "ssite_pos": str
            comma-separated list of cleavage site positions (0-indexed) within splice site
        "csite_inds": str
            comma-separate list of indices for each cleavage site in 3cs.txt
        "num_scens": str
            number of splice scenarios associated with each splice site site
        "ssite_start": str
            0-indexed start position of splice site
             * plus (+) strand: ssite_start is the first base of the splice site motif
             * minus (-) strand: ssite_end is the last base of the upsream intron
        "ssite_end": str
            0-indexed end position of splice site
             * plus (+) strand: ssite_end is the first base of the downstream exon
             * minus (-) strand: ssite_end the last base of the splice site motif
        "ssite_type": str
            splice site type, with information about number of trinucleotides and whether each is canonical (NAG)
        "ssite_seq": str
            chromosomal sequence of splice site
        "ssite_iflank_3ss": str
            100 bases upstream of splice site
        "ssite_eflank_3ss": str
            100 bases downstream of splice site
        "has_nag": int
            whether splice site contains a canonical NAG
             * 0: none of the trinucleotides in the splice site are a canonical NAG
             * 1: at least one trinucleotide in the splice site is a canonical NAG
    '''

    canon_chroms = [f"chr{c}" for c in list(range(1,23)) + ["X", "Y", "M"]]

    csites = {(c,s) : {} for c in canon_chroms for s in ["+", "-"]}
    ssites = {(c,s) : [] for c in canon_chroms for s in ["+", "-"]}

    # load cleavage sites
    cleavage_sites = pd.read_csv("data/3cs.txt", sep="\t")

    # add cleavage site positions + row index to correct cleavage site dictionaries
    for _,cs in cleavage_sites.iterrows():
        csites[(cs["chrom"], cs["strand"])][cs["down_start"]] = cs["index"]

    nags = ["AAG", "CAG", "GAG", "TAG"]

    # for each chromosome and strand:
    for chr_key in csites.keys():

        positions = list(csites[chr_key].keys()).copy()
        
        # for each csite position:
        for pos in csites[chr_key].keys():

            # if not already in a group:
            if pos in positions:
                group = {pos : csites[chr_key][pos]}
                positions.remove(pos)

                if True:
                    # search up
                    plus = 3
                    while ((pos + plus) in positions) or (get_3cs_tri(chr_key[0], pos+plus, chr_key[1], True) in nags):
                        # adjacent cleavage site
                        if (pos + plus) in positions:
                            group[pos+plus] = csites[chr_key][pos+plus]
                            positions.remove(pos+plus)
                        
                        # adjacent NAG
                        elif get_3cs_tri(chr_key[0], pos+plus, chr_key[1], True) in nags:
                            group[pos+plus] = -1
                        
                        plus += 3
                    
                    # search down
                    minus = 3
                    while ((pos - minus) in positions) or (get_3cs_tri(chr_key[0], pos-minus, chr_key[1], True) in nags):
                        # adjacent cleavage site
                        if (pos - minus) in positions:
                            group[pos-minus] = csites[chr_key][pos-minus]
                            positions.remove(pos-minus)
                        
                        # adjacent NAG
                        elif get_3cs_tri(chr_key[0], pos-minus, chr_key[1], True) in nags:
                            group[pos-minus] = -1

                        minus += 3
                else:
                    # search up
                    plus = 3
                    while (pos + plus) in positions:
                        group[pos+plus] = csites[chr_key][pos+plus]
                        positions.remove(pos+plus)
                        plus += 3

                    while get_3cs_tri(chr_key[0], pos+plus, chr_key[1], True) in nags:
                        group[pos+plus] = -1
                        pos += 3

                    # search down
                    minus = 3
                    while (pos - minus) in positions:
                        group[pos-minus] = csites[chr_key][pos-minus]
                        positions.remove(pos-minus)
                        minus += 3
                    
                    while get_3cs_tri(chr_key[0], pos-minus, chr_key[1], True) in nags:
                        group[pos-minus] = -1
                        minus += 3

                ssites[chr_key].append(group)

    # parse splice sites
    chrom = []
    strand = []

    num_tri = []
    num_csites = []
    csite_starts = []
    csite_pos = []
    csite_inds = []

    num_scens = []

    ssite_start = []
    ssite_end = []
    ssite_type = []

    ssite_seq = []
    iflank_3ss = []
    eflank_3ss = []
    has_nag = []

    # for each chromosome and strand:
    for chr_key,groups in ssites.items():

        # for each group:
        for group in groups:

            chrom.append(chr_key[0])
            strand.append(chr_key[1])

            # sort positions
            if chr_key[1] == "+":
                sorted_group = dict(sorted(group.items()))
            elif chr_key[1] == "-":
                sorted_group = dict(sorted(group.items(), reverse=True))

            cs = {s : i for s,i in sorted_group.items() if i != -1}

            num_tri.append(len(sorted_group))
            num_csites.append(len(cs))

            csite_starts.append(",".join([str(s) for s in cs.keys()]))
            csite_pos.append(",".join([str(list(sorted_group.keys()).index(s)) for s in cs.keys()]))
            csite_inds.append(",".join([str(i) for i in cs.values()]))
            num_scens.append(sum(cleavage_sites.iloc[list(cs.values())]["num_scens"]))

            if chr_key[1] == "+":
                start = list(sorted_group.keys())[0]-3
                end = list(sorted_group.keys())[-1]
            elif chr_key[1] == "-":
                start = list(sorted_group.keys())[0]+3
                end = list(sorted_group.keys())[-1]

            ssite_start.append(start)
            ssite_end.append(end)

            seqs = [get_3cs_tri(chr_key[0], p, chr_key[1]) for p in sorted_group.keys()]
            sequence = "".join(seqs)

            ssite_seq.append(sequence)
            ssite_type.append(get_sstype(sequence))

            iflank_3ss.append(get_ss_seq(chr_key[0], start, chr_key[1], 100, 0))
            eflank_3ss.append(get_ss_seq(chr_key[0], end, chr_key[1], 0, 100))

            has_nag.append(1 if sum([1 if s.upper() in ["AAG", "CAG", "GAG", "TAG"] else 0 for s in seqs]) > 0 else 0)

    # create DataFrame    
    ssites = pd.DataFrame({
        "chrom" : chrom,
        "strand" : strand,
        "num_tri" : num_tri,
        "num_csites" : num_csites,
        "csite_starts" : csite_starts,
        "csite_pos" : csite_pos,
        "csite_inds" : csite_inds,
        "num_scens" : num_scens,
        "ssite_start" : ssite_start,
        "ssite_end" : ssite_end,
        "ssite_type" : ssite_type,
        "ssite_seq" : ssite_seq,
        "ssite_iflank_3ss" : iflank_3ss,
        "ssite_eflank_3ss" : eflank_3ss,
        "has_nag" : has_nag
    })

    # add "index" column
    ssites = ssites.rename_axis("index").reset_index()

    # save DataFrame
    ssites.to_csv("data/3ss.txt", sep="\t", index=False)

# Identify uppercase
def identify_uppercase():

    ''' Identify splice sites whose entire motif is uppercase (does not fall within a repeat region)
    
    file: data/3ss.txt
    
        "uppercase": int
            whether any part of the splice site is part of a repetitive region
             * 0: one or more bases in the splice site motif are lowercase
             * 1: the entire splice site motif is uppercase    
    '''

    # load splice sites
    ssites = pd.read_csv("data/3ss.txt", sep="\t")

    uppercase = [0] * len(ssites)

    # for each splice site:
    for i,ss in ssites.iterrows():

        # mark if uppercase (no repeat masker)
        if ss["ssite_seq"] == ss["ssite_seq"].upper():
            uppercase[i] = 1

    # add column to DataFrame
    ssites["uppercase"] = uppercase

    # save DataFrame
    ssites.to_csv("data/3ss.txt", sep="\t", index=False)

# Extract NAGNAGs
def extract_nagnags():

    ''' Extract NAGNAG 3' cleavage & splice sites

    Save uppercase NAGNAG splice sites
        dst: data/nagnag_3ss.txt
    
    Save uppercase NAGNAG cleavage sites
        dst: data/nagnag_3cs.txt

        "ssite_start": str
            0-indexed start position, INCLUSIVE and irrespective of strand, of splice site
        "ssite_end": str
            0-indexed end position, EXCLUSIVE and irrespective of strand, of splice site
        "ssite_ind": int
            index of splice site in 3ss.txt
        "ssite_seq": str
            chromosomal sequence of splice site
        "ssite_down_seq": int
            first 3 bases after splice site
        "nagnag_stype": str {"PS", "DS", "AS"}
            cleavage type of NAGNAG
             * "PS": constitutive proximal
             * "DS": constitutive distal
             * "AS": alternative
        "nagnag_ctype": str {"CSP", "CSD", "ASP", "ASD"}
            cleavage type of NAGNAG
             * "CSP": constitutive proximal
             * "CSD": constitutive distal
             * "ASP": alternative proximal
             * "ASD": alternative distal
    '''

    # load splice sites
    ssites = pd.read_csv("data/3ss.txt", sep="\t")

    # load cleavage sites
    csites = pd.read_csv("data/3cs.txt", sep="\t")

    # isolate & save NAGNAG splice sites
    nagnag_ss = ssites[(ssites["uppercase"] == 1) & (ssites["ssite_type"] == "2C")]
    nagnag_ss.to_csv("data/nagnag_3ss.txt", sep="\t", index=False)

    nagnag_ss.reset_index(inplace=True)

    ssite_start = [0] * len(csites)
    ssite_end = [0] * len(csites)
    ssite_ind = [0] * len(csites)
    ssite_seq = [""] * len(csites)
    ssite_down_seq = [""] * len(csites)
    nagnag_stype = [""] * len(csites)
    nagnag_ctype = [""] * len(csites)

    # for each NAGNAG splice site:
    for _,nagnag in nagnag_ss.iterrows():

        # for each cleavage site:
        for p,j in zip(nagnag["csite_pos"].split(","), nagnag["csite_inds"].split(",")):
            ssite_seq[int(j)] = nagnag["ssite_seq"]
            ssite_start[int(j)] = nagnag["ssite_start"]
            ssite_end[int(j)] = nagnag["ssite_end"]
            ssite_ind[int(j)] = nagnag["index"]
            ssite_down_seq[int(j)] = str(nagnag["ssite_eflank_3ss"][:3])

            # alternatively-spliced NAGNAG
            if nagnag["csite_pos"] == "0,1":
                nagnag_stype[int(j)] = "AS"

                # proximal splice site
                if p == "0":
                    nagnag_ctype[int(j)] = "ASP"

                # distal cleavage site
                elif p == "1":
                    nagnag_ctype[int(j)] = "ASD"

            # constitutive proximally-spliced NAGNAG
            elif nagnag["csite_pos"] == "0":
                nagnag_stype[int(j)] = "PS"
                nagnag_ctype[int(j)] = "CSP"
            
            # constitutive distally-spliced NAGNAG
            elif nagnag["csite_pos"] == "1":
                nagnag_stype[int(j)] = "DS"
                nagnag_ctype[int(j)] = "CSD"

    # add columns to cleavage site DataFrame
    csites["ssite_start"] = ssite_start
    csites["ssite_end"] = ssite_end
    csites["ssite_ind"] = ssite_ind
    csites["ssite_seq"] = ssite_seq
    csites["ssite_down_seq"] = ssite_down_seq
    csites["nagnag_stype"] = nagnag_stype
    csites["nagnag_ctype"] = nagnag_ctype

    # get cleavage site row indices
    csite_inds = [int(ind) for i in range(len(nagnag_ss)) for ind in nagnag_ss["csite_inds"][i].split(",")]

    # isolate & save NAGNAG cleavage sites
    nagnag_cs = csites.iloc[csite_inds]
    nagnag_cs.to_csv("data/nagnag_3cs.txt", sep="\t", index=False)

# Extract NAGNAG splice scenarios
def split_nagnag_scenarios():

    ''' Split NAGNAGs into splice scenarios
    
    src: data/nagnag_3ss.txt; data/nagnag_3cs.txt
    
    dst: data/nagnag_scens.txt
    
        "index": int
            0-indexed row number in DataFrame
        "chrom": str
            chromosome on which splice scenario is found
            format: chr{1, 2, ... , 21, 22, X, Y, M}
        "down_start": int
            0-indexed chromosome position of start of downstream exon
             * plus strand: first base of downstream exon 
             * minus strand: last base of upstream intron                
        "strand": str ("+", "-")
            strand of transcript in which splice site is found
        "rtype": str {"CDS", "5'UTR", "3'UTR", "ncRNA"}
            comma-separated list, 1 for each splice scenario, of the RNA type of cleavage site and transcript
             * "CDS: found in a coding (NM, XM) transcript, within the coding sequence
             * "5'UTR": found in a coding (NM, XM) transcript, upstream of the coding sequence
             * "3'UTR": found in a coding (NM, XM) transcript, downstream of the coding sequence
             * "ncRNA":  found in an non-coding (NR, XR) transcript
        "up_end": int
            0-indexed chromosome position of end of upstream exon
             * plus strand: first base of downstream intron 
             * minus strand: last base of upstream exon
        "frame": int {-1, 0, 1, 2}
            the frame/phase of the downstream exon
             * -1: downstream exon is in a UTR or in a non-coding RNA transcript
             * 0: downstream exon is the first base of its codon
             * 1: downstream exon is the second base of its codon
             * 2: downstream exon is the third base of its codon
        "accession": str
            comma-separated list of forward-slash-separated lists of accession numbers associated with this splice scenario
        "num_isoforms": int
            comma-separated list, 1 for each splice scenario, of the number of transcripts associated with each splice scenario
        "scen_up_seq": str
            comma-separate list, 1 for each splice scenario, of the last 3 bases of the upstream exon
        "ssite_seq": str
            chromosomal sequence of splice site
        "ssite_start": str
            0-indexed start position, INCLUSIVE and irrespective of strand, of splice site
        "ssite_end": str
            0-indexed end position, EXCLUSIVE and irrespective of strand, of splice site
        "ssite_down_seq": int
            first 3 bases after splice site
        "nagnag_stype": str {"PS", "DS", "AS"}
            cleavage type of NAGNAG
             * "PS": constitutive proximal
             * "DS": constitutive distal
             * "AS": alternative
        "nagnag_ctype": str {"CSP", "CSD", "ASP", "ASD"}
            cleavage type of NAGNAG scenario
             * "CSP": constitutive proximal
             * "CSD": constitutive distal
             * "ASP": alternative proximal
             * "ASD": alternative distal
        "csite_ind": int
            index of associated cleavage site in nagnag_3cs.txt
        "ssite_ind": int
            index of associated cleavage site in nagnag_3ss.txt

    Add scen index to NAGNAG splice sites
        file: data/nagnag_3ss.txt

        "scen_ind": str
            comma-separated list, 1 for each cleavage site, of forward-slash-separated lists of indices, one for each splice scenario, in nagnag_scens.txt

    Add scen index to NAGNAG cleavage sites
        file: data/nagnag_3cs.txt

        "scen_ind": str
            comma-separated list of indices of splice scenarios in nagnag_scens.txt
    '''

    # load splice sites
    ssites = pd.read_csv("data/nagnag_3ss.txt", sep="\t", index_col=0)

    scen_inds_ss = [""] * len(ssites)

    # load cleavage sites
    csites = pd.read_csv("data/nagnag_3cs.txt", sep="\t", index_col=0)

    scen_inds_cs = [""] * len(csites)
    
    chrom = [""] * sum(csites["num_scens"])
    down_start = [0] * sum(csites["num_scens"])
    strand = [""] * sum(csites["num_scens"])
    rtype = [""] * sum(csites["num_scens"])
    up_end = [0] * sum(csites["num_scens"])
    frame = [0] * sum(csites["num_scens"])
    accession = [""] * sum(csites["num_scens"])
    n_isoforms = [0] * sum(csites["num_scens"])
    up_seq = [""] * sum(csites["num_scens"])
    ssite_seq = [""] * sum(csites["num_scens"])
    ssite_start = [0] * sum(csites["num_scens"])
    ssite_end = [0] * sum(csites["num_scens"])
    ssite_down_seq = [""] * sum(csites["num_scens"])
    nagnag_stype = [""] * sum(csites["num_scens"])
    nagnag_ctype = [""] * sum(csites["num_scens"])
    csite_ind = [0] * sum(csites["num_scens"])
    ssite_ind = [0] * sum(csites["num_scens"])

    s = 0

    # for each cleavage site:
    for i,(c,cs) in enumerate(csites.iterrows()):

        scen_inds = []

        cs_rt = cs["rtype"].split(",")
        cs_up_end = cs["up_end"].split(",")
        cs_fr = cs["frame"].split(",")
        cs_acc = cs["accession"].split(",")
        cs_n_iso = cs["n_isoforms"].split(",")
        cs_up_seq = cs["scen_eflank_5ss"].split(",")

        for rt,ue,f,a,n,us in zip(cs_rt, cs_up_end, cs_fr, cs_acc, cs_n_iso, cs_up_seq):
            chrom[s] = cs["chrom"]
            down_start[s] = cs["down_start"]
            strand[s] = cs["strand"]
            rtype[s] = rt
            up_end[s] = ue
            frame[s] = f
            accession[s] = a
            n_isoforms[s] = n
            up_seq[s] = us[-3:]
            ssite_seq[s] = cs["ssite_seq"]
            ssite_start[s] = cs["ssite_start"]
            ssite_end[s] = cs["ssite_end"]
            ssite_down_seq[s] = cs["ssite_down_seq"]
            nagnag_stype[s] = cs["nagnag_stype"]
            nagnag_ctype[s] = cs["nagnag_ctype"]
            csite_ind[s] = c
            ssite_ind[s] = cs["ssite_ind"]
            
            scen_inds.append(s)

            s += 1
        
        scen_inds_cs[i] = ",".join([str(si) for si in scen_inds])
    
    # create DataFrame    
    scens = pd.DataFrame({
        "chrom" : chrom,
        "down_start" : down_start,
        "strand" : strand,
        "rtype" : rtype,
        "up_end" : up_end,
        "frame" : frame,
        "accession" : accession,
        "n_isoforms" : n_isoforms,
        "scen_up_seq" : up_seq,
        "ssite_seq" : ssite_seq,
        "ssite_start" : ssite_start,
        "ssite_end" : ssite_end,
        "ssite_down_seq" : ssite_down_seq,
        "nagnag_stype" : nagnag_stype,
        "nagnag_ctype" : nagnag_ctype,
        "csite_ind" : csite_ind,
        "ssite_ind" : ssite_ind
    })

    # add "index" column
    scens = scens.rename_axis("index").reset_index()

    # save DataFrame
    scens.to_csv("data/nagnag_scens.txt", sep="\t", index=False)

    # add "scen_inds" to cleavage sites
    csites["scen_inds"] = scen_inds_cs
    csites.to_csv("data/nagnag_3cs.txt", sep="\t", index=True)

    # for each splice site:
    for i,(_,ss) in enumerate(ssites.iterrows()):
        cs_inds = [int(ci) for ci in ss["csite_inds"].split(",")]
        scen_inds_ss[i] = ",".join(["/".join(csites["scen_inds"][ci].split(",")) for ci in cs_inds])

    # add "scen_inds" to splice sites
    ssites["scen_inds"] = scen_inds_ss
    ssites.to_csv("data/nagnag_3ss.txt", sep="\t", index=True)

# Add NAGNAG splice scenario flanking sequences
def add_scen_flanks():

    ''' Get sequences flanking NAGNAG splice scenarios

    file: data/nagnag_scens.txt
    
        "eflank_5ss": str
            100 exonic bases upstream of the 5' cleavage site
        "iflank_5ss": str
            100 intronic bases downstream of the 5' cleavage site
        "iflank_3ss": str
            100 intronic bases upstream of the NAGNAG motif
        "eflank_3ss": str
            100 exonic bases downstream of the NAGNAG motif
    '''

    # load splice sites
    scens = pd.read_csv("data/nagnag_scens.txt", sep="\t", index_col=0)

    iflank_3ss = [""] * len(scens)
    eflank_3ss = [""] * len(scens)

    iflank_5ss = [""] * len(scens)
    eflank_5ss = [""] * len(scens)

    # for each splice site
    for i,(_,scen) in enumerate(scens.iterrows()):
        iflank_3ss[i] = get_ss_seq(scen["chrom"], scen["ssite_start"], scen["strand"], 100, 0)
        eflank_3ss[i] = get_ss_seq(scen["chrom"], scen["ssite_end"], scen["strand"], 0, 100)

        iflank_5ss[i] = get_ss_seq(scen["chrom"], scen["up_end"], scen["strand"], 0, 100)
        eflank_5ss[i] = get_ss_seq(scen["chrom"], scen["up_end"], scen["strand"], 100, 0)

    # add columns to DataFrame
    scens["scen_eflank_5ss"] = eflank_5ss
    scens["scen_iflank_5ss"] = iflank_5ss
    scens["ssite_iflank_3ss"] = iflank_3ss
    scens["ssite_eflank_3ss"] = eflank_3ss

    # save DataFrame
    scens.to_csv("data/nagnag_scens.txt", sep="\t", index=True)
