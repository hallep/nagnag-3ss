''' Parse NCBI RefSeq transcript annotation data 
    and extract splice sites, cleavage sites, and splice scenarios '''

from lib import ROOT, pd
from sequence import get_ss_seq, get_cs_seq, is_nag, get_sstype

# Get cleavage sites
def get_csites():

    ''' Identify all unique 3' cleavage sites
    
    **Source:** `src/ncbiRefSeq.tx`    
    
    **Destination:** `sites/3cs.txt`
    * index (int): 0-indexed row number in DataFrame
    * chrom (str): chromosome on which 3' cleavage site is found, formatted as 'chr#'
    * down_start (int): 0-indexed chromosome position of start of downstream exon
        * plus strand: first base of downstream exon 
        * minus strand: last base of upstream intron
    * strand (str): strand of transcript in which cleavage site is found
        * possible values: "+", "-"
    * num_scens (int): number of splice scenarios associated with each cleavage site
    * up_end (str): comma-separate list, 1 for each splice scenario, of the 0-indexed chromosome position of end of upstream exon
        * plus strand: first base of downstream intron 
        * minus strand: last base of upstream exon
    * rtype (str): comma-separated list, 1 for each splice scenario, of the RNA type of cleavage site and transcript
        * "CDS: found in a coding (NM, XM) transcript, within the coding sequence
        * "5'UTR": found in a coding (NM, XM) transcript, upstream of the coding sequence
        * "3'UTR": found in a coding (NM, XM) transcript, downstream of the coding sequence
        * "ncRNA":  found in an non-coding (NR, XR) transcript
    * phase (str): comma-separated list, 1 for each splice scenario, of the phase (frame) of downstream exon
        * -1: downstream exon is in a UTR or in a non-coding RNA transcript
        * 0: downstream exon is the first base of its codon
        * 1: downstream exon is the second base of its codon
        * 2: downstream exon is the third base of its codon
    * accession (str): comma-separated list, 1 for each splice scenario, of forward-slash-separated lists of accession numbers
    * n_isoforms (str): comma-separated list, 1 for each splice scenario, of the number of transcripts associated with each splice scenario
    * csite_seq (str): last 3 bases of intron
    * scen_eflank_5ss (str): comma-separated list, 1 for each splice scenario, of the last 100 bases of the upstream exon
    * scen_iflank_5ss (str): comma-separated list, 1 for each splice scenario, of the first 100 bases of the intron following the upstream exon
    * csite_iflank_3ss (str): last 100 bases of the intron (including csite_seq)
    * csite_eflank_3ss (str): first 100 bases of the downstream exon
    '''

    # load transcript data table
    columns = ["bin", "name", "chrom", "strand", "txStart", "txEnd",
               "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
               "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"]
    ncbi_table = pd.read_csv(f"{ROOT}/src/ncbiRefSeq.txt", sep="\t", header=None, names=columns)

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
    d = {
        "chrom" : [""] * len(cs3),
        "down_start" : [""] * len(cs3),
        "strand" : [""] * len(cs3),
        "num_scens" : [""] * len(cs3),
        "rtype" : [""] * len(cs3),
        "up_end" : [""] * len(cs3),
        "phase" : [""] * len(cs3),
        "accession" : [""] * len(cs3),
        "n_isoforms" : [""] * len(cs3),
        "csite_seq" : [""] * len(cs3),
        "scen_eflank_5ss" : [""] * len(cs3),
        "scen_iflank_5ss" : [""] * len(cs3),
        "csite_iflank_3ss" : [""] * len(cs3),
        "csite_eflank_3ss" : [""] * len(cs3),
    }

    # for each cleavage site:
    for i,(csite,scens) in enumerate(cs3.items()):

        # coordinates
        d["chrom"][i] = csite[0]
        d["down_start"][i] = csite[1]
        d["strand"][i] = csite[2]

        # scenarios
        d["num_scens"][i] = len(scens)

        r,f,u = zip(*list(scens))

        d["rtype"][i] = ",".join(r)
        d["phase"][i] = ",".join(f)
        d["up_end"][i] = ",".join(map(str, u))

        d["accession"][i] = ",".join(["/".join(acc_nums) for acc_nums in scens.values()])
        d["n_isoforms"][i] = ",".join([str(len(acc_nums)) for acc_nums in scens.values()])
        
        # 3' splice site sequences
        d["csite_seq"][i] = get_cs_seq(csite[0], csite[1], csite[2])
        d["csite_iflank_3ss"][i] = get_ss_seq(csite[0], csite[1], csite[2], 100, 0)
        d["csite_eflank_3ss"][i] = get_ss_seq(csite[0], csite[1], csite[2], 0, 100)

        # 5' splice site sequences
        fe = [""] * len(scens)
        fi = [""] * len(scens)

        for j,s in enumerate(u):
            fe[j] = get_ss_seq(csite[0], s, csite[2], 100, 0)
            fi[j] = get_ss_seq(csite[0], s, csite[2], 0, 100)

        d["scen_eflank_5ss"][i] = ",".join(fe)
        d["scen_iflank_5ss"][i] = ",".join(fi)

    # DataFrame    
    csites = pd.DataFrame(d)
    csites.to_csv(f"{ROOT}/sites/3cs.txt", sep="\t", index_label="index")

# Group cleavage sites
def group_csites():

    ''' Group 3' cleavage sites into 3' splice sites
    
    **Source:** `sites/3cs.txt`
    
    **Destination:** `sites/3ss.txt`
    * index (int): 0-indexed row number in DataFrame
    * chrom (str): chromosome on which 3' splice site is found, formatted as 'chr#'
    * strand (str): strand of transcript in which cleavage site is found
        * possible values: "+", "-"
    * num_tri (int): number of cleavage sites + adjacent NAGs in the splice site
    * num_csites (int): number of cleavage sites in the splice site
    * csite_starts (str): comma-separated list of cleavage site down_start positions
    * csite_pos (str): comma-separated list of cleavage site positions (0-indexed) within splice site
    * csite_inds (str): comma-separate list of indices for each cleavage site in 3cs.txt
    * num_scens (str): number of splice scenarios associated with each splice site
    * ssite_start (int): 0-indexed start position of splice site (inclusive)
        * plus (+) strand: first base of the splice site motif
        * minus (-) strand: last base of the upstream intron
    * ssite_end (int): 0-indexed end position of splice site (exclusive)
        * plus (+) strand: first base of the downstream exon
        * minus (-) strand: last base of the splice site motif
    * ssite_type (str): splice site type, with information about number of trinucleotides and whether each is canonical (NAG)
    * ssite_seq (str): genomic sequence of splice site
    * ssite_iflank_3ss (str): 100 bases upstream of (and not including) splice site
    * ssite_eflank_3ss (str): 100 bases downstream of (and not including) splice site
    * has_nag (int): whether splice site contains a canonical NAG
        * 0: none of the trinucleotides in the splice site are a canonical NAG
        * 1: at least one trinucleotide in the splice site is a canonical NAG
    * uppercase (int): whether the entire splice site sequence is uppercase

    **Destination:** `sites/3cs.txt`
    * ssite_ind (int): index of splice site in 3ss.txt
    '''

    canon_chroms = [f"chr{c}" for c in list(range(1,23)) + ["X", "Y", "M"]]

    csites = {(c,s) : {} for c in canon_chroms for s in ["+", "-"]}
    ssites = {(c,s) : [] for c in canon_chroms for s in ["+", "-"]}

    # load cleavage sites
    cleavage_sites = pd.read_csv(f"{ROOT}/sites/3cs.txt", sep="\t")

    # add cleavage site positions + row index to correct cleavage site dictionaries
    for _,cs in cleavage_sites.iterrows():
        csites[(cs["chrom"], cs["strand"])][cs["down_start"]] = cs["index"]

    # for each chromosome and strand:
    for chr_key in csites.keys():

        unused = set(csites[chr_key])
                
        # for each cleavage site position:
        for pos in csites[chr_key].keys():

            # if not already in a group:
            if pos in unused:
                group = {pos : csites[chr_key][pos]}
                unused.remove(pos)

                # search up
                plus = 3
                while ((pos + plus) in unused) or is_nag(get_cs_seq(chr_key[0], pos+plus, chr_key[1])):
                    
                    # adjacent cleavage site
                    if (pos + plus) in unused:
                        group[pos+plus] = csites[chr_key][pos+plus]
                        unused.remove(pos+plus)
                    
                    # adjacent NAG
                    elif is_nag(get_cs_seq(chr_key[0], pos+plus, chr_key[1])):
                        group[pos+plus] = -1
                    
                    plus += 3
                
                # search down
                minus = 3
                while ((pos - minus) in unused) or is_nag(get_cs_seq(chr_key[0], pos-minus, chr_key[1])):
                    
                    # adjacent cleavage site
                    if (pos - minus) in unused:
                        group[pos-minus] = csites[chr_key][pos-minus]
                        unused.remove(pos-minus)
                    
                    # adjacent NAG
                    elif is_nag(get_cs_seq(chr_key[0], pos-minus, chr_key[1])):
                        group[pos-minus] = -1

                    minus += 3

                # save group
                ssites[chr_key].append(group)

    # parse splice sites
    d = {
        "chrom" : [],
        "strand" : [],
        "num_tri" : [],
        "num_csites" : [],
        "csite_starts" : [],
        "csite_pos" : [],
        "csite_inds" : [],
        "num_scens" : [],
        "ssite_start" : [],
        "ssite_end" : [],
        "ssite_type" : [],
        "ssite_seq" : [],
        "ssite_iflank_3ss" : [],
        "ssite_eflank_3ss" : [],
        "has_nag" : [],
        "uppercase" : [],
    }

    # for each chromosome and strand:
    for chr_key,groups in ssites.items():

        # for each group:
        for group in groups:

            d["chrom"].append(chr_key[0])
            d["strand"].append(chr_key[1])

            # sort positions
            if chr_key[1] == "+":
                sorted_group = dict(sorted(group.items()))
            elif chr_key[1] == "-":
                sorted_group = dict(sorted(group.items(), reverse=True))

            cs = {s : i for s,i in sorted_group.items() if i != -1}

            d["num_tri"].append(len(sorted_group))
            d["num_csites"].append(len(cs))

            d["csite_starts"].append(",".join([str(s) for s in cs.keys()]))
            d["csite_pos"].append(",".join([str(list(sorted_group.keys()).index(s)) for s in cs.keys()]))
            d["csite_inds"].append(",".join([str(i) for i in cs.values()]))
            d["num_scens"].append(cleavage_sites.iloc[list(cs.values())]["num_scens"].sum())

            if chr_key[1] == "+":
                start = list(sorted_group.keys())[0]-3
                end = list(sorted_group.keys())[-1]
            elif chr_key[1] == "-":
                start = list(sorted_group.keys())[0]+3
                end = list(sorted_group.keys())[-1]

            d["ssite_start"].append(start)
            d["ssite_end"].append(end)

            seqs = [get_cs_seq(chr_key[0], p, chr_key[1]) for p in sorted_group.keys()]
            sequence = "".join(seqs)

            d["ssite_type"].append(get_sstype(sequence))

            d["ssite_seq"].append(sequence)
            d["ssite_iflank_3ss"].append(get_ss_seq(chr_key[0], start, chr_key[1], 100, 0))
            d["ssite_eflank_3ss"].append(get_ss_seq(chr_key[0], end, chr_key[1], 0, 100))

            d["has_nag"].append(int(any(is_nag(s) for s in seqs)))
            d["uppercase"].append(int(sequence.isupper()))

    # splice site DataFrame
    ssites = pd.DataFrame(d)
    ssites.to_csv(f"{ROOT}/sites/3ss.txt", sep="\t", index_label="index")

    # splice site indices
    ss_inds = [0] * len(cleavage_sites)

    for si,cis in ssites["csite_inds"].items():
        for ci in map(int, cis.split(",")):
            ss_inds[ci] = si

    cleavage_sites["ssite_ind"] = ss_inds
    cleavage_sites.to_csv(f"{ROOT}/sites/3cs.txt", sep="\t", index=False)

# Isolate NAGNAGs
def isolate_nagnags():

    ''' Isolate (uppercase) NAGNAG 3' cleavage & splice sites

    **Source:**
    * splice sites: `sites/3ss.txt`
    * cleavage sites: `sites/3cs.txt`

    **Destination:** `sites/nagnag_3ss.txt`

    **Destination:** `sites/nagnag_3cs.txt`
    * ssite_start (int): 0-indexed start position, INCLUSIVE and irrespective of strand, of splice site
    * ssite_end (int): 0-indexed end position, EXCLUSIVE and irrespective of strand, of splice site
    * ssite_seq (str): genomic sequence of splice site
    * ssite_down_seq (int): first 3 bases after splice site
    * ssite_stype (str): splice type of NAGNAG splice site
        * "PS": constitutive proximal
        * "DS": constitutive distal
        * "AS": alternative
    * csite_stype (str): splice type of cleavage site
        * "CSP": constitutive proximal
        * "CSD": constitutive distal
        * "ASP": alternative proximal
        * "ASD": alternative distal
    '''

    # load sites
    ssites = pd.read_csv(f"{ROOT}/sites/3ss.txt", sep="\t")
    csites = pd.read_csv(f"{ROOT}/sites/3cs.txt", sep="\t")

    # isolate NAGNAG splice sites
    nagnag_ss = ssites[(ssites["uppercase"] == 1) & (ssites["ssite_type"] == "2C")]
    nagnag_ss.to_csv(f"{ROOT}/sites/nagnag_3ss.txt", sep="\t", index=False)

    # isolate cleavage sites
    nagnag_ss.reset_index(inplace=True)

    ssite_start = [0] * len(csites)
    ssite_end = [0] * len(csites)
    ssite_ind = [0] * len(csites)
    ssite_seq = [""] * len(csites)
    ssite_down_seq = [""] * len(csites)
    ss_stype = [""] * len(csites)
    cs_stype = [""] * len(csites)

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
                ss_stype[int(j)] = "AS"

                # proximal splice site
                if p == "0":
                    cs_stype[int(j)] = "ASP"

                # distal cleavage site
                elif p == "1":
                    cs_stype[int(j)] = "ASD"

            # constitutive proximally-spliced NAGNAG
            elif nagnag["csite_pos"] == "0":
                ss_stype[int(j)] = "PS"
                cs_stype[int(j)] = "CSP"
            
            # constitutive distally-spliced NAGNAG
            elif nagnag["csite_pos"] == "1":
                ss_stype[int(j)] = "DS"
                cs_stype[int(j)] = "CSD"

    # add columns to cleavage site DataFrame
    csites["ssite_start"] = ssite_start
    csites["ssite_end"] = ssite_end
    csites["ssite_ind"] = ssite_ind
    csites["ssite_seq"] = ssite_seq
    csites["ssite_down_seq"] = ssite_down_seq
    csites["ssite_stype"] = ss_stype
    csites["csite_stype"] = cs_stype

    # isolate & save NAGNAG cleavage sites
    csite_inds = [int(i) for inds in nagnag_ss["csite_inds"] for i in inds.split(",")]
    nagnag_cs = csites.iloc[csite_inds]
    nagnag_cs.to_csv(f"{ROOT}/sites/nagnag_3cs.txt", sep="\t", index=False)

# Parse NAGNAG splice scenarios
def parse_nagnag_scenarios():

    ''' Split NAGNAGs into splice scenarios

    **Source:**
    * splice sites: `sites/nagnag_3ss.txt`
    * cleavage sites: `sites/nagnag_3cs.txt`

    **Destination:** `sites/nagnag_scens.txt`
    * index (int): 0-indexed row number in DataFrame
    * chrom (str): chromosome on which splice scenario is found, formatted as 'chr#'
    * ssite_start (int): 0-indexed start position, INCLUSIVE and irrespective of strand, of splice site
    * ssite_end (int): 0-indexed end position, EXCLUSIVE and irrespective of strand, of splice site
    * strand (str): strand of transcript in which splice scenario is found
        * possible values: "+", "-"
    * up_end (int): 0-indexed chromosome position of end of upstream exon
        * plus strand: first base of downstream intron
        * minus strand: last base of upstream exon
    * rtype (str): RNA type of splice scenario
        * "CDS: found in a coding (NM, XM) transcript, within the coding sequence
        * "5'UTR": found in a coding (NM, XM) transcript, upstream of the coding sequence
        * "3'UTR": found in a coding (NM, XM) transcript, downstream of the coding sequence
        * "ncRNA":  found in an non-coding (NR, XR) transcript
    * phase (int): phase (frame) of downstream exon
        * -1: downstream exon is in a UTR or in a non-coding RNA transcript
        * 0: downstream exon is the first base of its codon
        * 1: downstream exon is the second base of its codon
        * 2: downstream exon is the third base of its codon
    * accession (str): comma-separated list, one for each cleavage site, of forward-slash-separated lists of accession numbers
    * n_isoforms (str): comma-separated list, one for each cleavage site, of number of transcripts associated with splice scenario
    * up_seq (str): last 3 bases of the upstream exon
    * ssite_seq (str): genomic sequence of splice site
    * down_seq (str): first 3 bases after splice site
    * splice_type (str): splice type of splice scenario
        * "PS": constitutive proximal
        * "DS": constitutive distal
        * "AS": alternative
    * csite_inds (str): comma-separated list of indices of associated cleavage sites in nagnag_3cs.txt
    * ssite_ind (int): index of associated cleavage site in nagnag_3ss.txt

    **Destination:** `sites/nagnag_3ss.txt`
    * scen_inds (str): comma-separated list of indices of associated splice scenarios in nagnag_scens.txt

    **Destination:** `sites/nagnag_3cs.txt`
    * scen_ind (str): comma-separated list, 1 for each scenario, of indices of associated splice scenarios in nagnag_scens.txt
    '''

    # load splice sites
    ssites = pd.read_csv(f"{ROOT}/sites/nagnag_3ss.txt", sep="\t", index_col=0)
    scen_inds_ss = {i : [] for i in ssites.index}

    # load cleavage sites
    csites = pd.read_csv(f"{ROOT}/sites/nagnag_3cs.txt", sep="\t", index_col=0)
    scen_inds_cs = {i : ["" for _ in range(n)] for i,n in csites["num_scens"].items()}
    
    d = {
        "chrom" : [],
        "ssite_start" : [],
        "ssite_end" : [],
        "strand" : [],

        "up_end" : [],
        "rtype" : [],
        "phase" : [],
        "accession" : [],
        "n_isoforms" : [],
        
        "up_seq" : [],
        "ssite_seq" : [],
        "down_seq" : [],
        
        "splice_type" : [],
        "csite_inds" : [],
        "ssite_ind" : [],
    }

    x = 0
    
    # for each cleavage site:
    for si,ssite in ssites.iterrows():
        
        # cleavage sites
        cs = {int(p) : int(c) for p,c in zip(ssite["csite_pos"].split(","), ssite["csite_inds"].split(","))}

        # splice instances: (up_end, rtype, phase, up_seq)
        insts = {0 : [], 1 : []}
        insts.update({cp : {(u,r,p,s) : (j,a,n) for j,(u,r,p,s,a,n) in enumerate(zip(csites.loc[ci]["up_end"].split(","),
                                                                                     csites.loc[ci]["rtype"].split(","),
                                                                                     csites.loc[ci]["phase"].split(","),
                                                                                     csites.loc[ci]["scen_eflank_5ss"].split(","),
                                                                                     csites.loc[ci]["accession"].split(","),
                                                                                     csites.loc[ci]["n_isoforms"].split(",")))}
                            for cp,ci in cs.items()})
    
        # for all splice scenarios:
        for scen in set(list(insts[0]) + list(insts[1])):
            
            # location information (first base after motif)
            d["chrom"].append(ssite["chrom"])
            d["ssite_start"].append(ssite["ssite_start"])
            d["ssite_end"].append(ssite["ssite_end"])
            d["strand"].append(ssite["strand"])
            
            # scenario information
            d["up_end"].append(scen[0])
            d["rtype"].append(scen[1])
            d["phase"].append(scen[2])

            # sequences
            d["up_seq"].append(scen[3][-3:])
            d["ssite_seq"].append(ssite["ssite_seq"])
            d["down_seq"].append(ssite["ssite_eflank_3ss"][:3])

            # alternatively-spliced
            if (scen in insts[0]) and (scen in insts[1]):
                d["splice_type"].append("AS")
                pos = [0, 1]
            
            # proximally-spliced
            elif (scen in insts[0]):
                d["splice_type"].append("PS")
                pos = [0]

            # distally-spliced
            elif (scen in insts[1]):
                d["splice_type"].append("DS")
                pos = [1]

            d["accession"].append(",".join([insts[p][scen][1] for p in pos]))
            d["n_isoforms"].append(",".join([insts[p][scen][2] for p in pos]))

            # indices
            d["csite_inds"].append(",".join(map(str, [cs[p] for p in pos])))
            d["ssite_ind"].append(si)

            scen_inds_ss[si].append(x)
            for p in pos:
                scen_inds_cs[cs[p]][insts[p][scen][0]] = x

            x += 1

    # DataFrame    
    scens = pd.DataFrame(d)
    scens.to_csv(f"{ROOT}/sites/nagnag_scens.txt", sep="\t", index_label="index")

    print(scens)

    # splice scenario indices
    ssites["scen_inds"] = [",".join(map(str, inds)) for inds in scen_inds_ss.values()]
    ssites.to_csv(f"{ROOT}/sites/nagnag_3ss.txt", sep="\t", index=True)

    csites["scen_inds"] = [",".join(map(str, inds)) for inds in scen_inds_cs.values()]
    csites.to_csv(f"{ROOT}/sites/nagnag_3cs.txt", sep="\t", index=True)

    print(ssites)
    print(csites)

get_csites()
group_csites()
isolate_nagnags()
parse_nagnag_scenarios()
