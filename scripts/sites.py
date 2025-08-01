''' Parse NCBI RefSeq transcript annotation data 
    and extract splice sites, cleavage sites, and splice scenarios '''

from lib import ROOT, re, pd
from sequence import ss_seq

# Get cleavage sites
def get_cleavage_sites():

    ''' Identify all unique 3' cleavage sites
    
    Source: src/ncbiRefSeq.tx    
    
    Destination: sites/3cs.txt
    * index (int): 0-indexed row number in DataFrame
    * chrom (str): chromosome on which 3' cleavage site is found, formatted as 'chr#'
    * down_start (int): 0-indexed chromosome position of start of downstream exon
        * plus strand: first base of downstream exon 
        * minus strand: last base of upstream intron
    * strand (str): strand of transcript in which cleavage site is found
        * possible values: "+", "-"
    * num_scens (int): number of splice scenarios associated with each cleavage site
    * rtype (str): comma-separated list, 1 for each splice scenario, of the RNA type of cleavage site and transcript
        * "CDS: found in a coding (NM, XM) transcript, within the coding sequence
        * "5'UTR": found in a coding (NM, XM) transcript, upstream of the coding sequence
        * "3'UTR": found in a coding (NM, XM) transcript, downstream of the coding sequence
        * "ncRNA":  found in an non-coding (NR, XR) transcript
    * up_end (str): comma-separate list, 1 for each splice scenario, of the 0-indexed chromosome position of end of upstream exon
        * plus strand: first base of downstream intron 
        * minus strand: last base of upstream exon
    * frame (str): comma-separated list, 1 for each splice scenario, of the frame/phase of downstream exon
        * -1: downstream exon is in a UTR or in a non-coding RNA transcript
        * 0: downstream exon is the first base of its codon
        * 1: downstream exon is the second base of its codon
        * 2: downstream exon is the third base of its codon
    * accession (str): comma-separated list, 1 for each splice scenario, of forward-slash-separated lists of accession numbers
    * n_isoforms (str): comma-separated list, 1 for each splice scenario, of the number of transcripts associated with each splice scenario
    * csite_seq (str): last 3 bases of intron
    * scen_eflank_5ss (str): comma-separated list, 1 for each splice scenario, of the last 100 bases of the upstream exon
    * scen_iflank_5ss (str): comma-separated list, 1 for each splice scenario, of the first 100 bases of the intron following the upstream exon
    * csite_iflank_3ss (str): last 100 bases of the intron (including csite_seqs)
    * csite_eflank_3ss (str): first 100 bases of the downstream exon
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
    d = {
        "chrom" : [""] * len(cs3),
        "down_start" : [""] * len(cs3),
        "strand" : [""] * len(cs3),
        "num_scens" : [""] * len(cs3),
        "rtype" : [""] * len(cs3),
        "up_end" : [""] * len(cs3),
        "frame" : [""] * len(cs3),
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
        d["frame"][i] = ",".join(f)
        d["up_end"][i] = ",".join(map(str, u))

        d["accession"][i] = ",".join(["/".join(acc_nums) for acc_nums in scens.values()])
        d["n_isoforms"][i] = ",".join([str(len(acc_nums)) for acc_nums in scens.values()])
        
        # 3' splice site sequences
        d["csite_seq"][i] = ss_seq(csite[0], csite[1], csite[2], 3, 0)
        d["csite_iflank_3ss"][i] = ss_seq(csite[0], csite[1], csite[2], 100, 0)
        d["csite_eflank_3ss"][i] = ss_seq(csite[0], csite[1], csite[2], 0, 100)

        # 5' splice site sequences
        fe = [""] * len(scens)
        fi = [""] * len(scens)

        for j,s in enumerate(u):
            fe[j] = ss_seq(csite[0], s, csite[2], 100, 0)
            fi[j] = ss_seq(csite[0], s, csite[2], 0, 100)

        d["scen_eflank_5ss"][i] = ",".join(fe)
        d["scen_iflank_5ss"][i] = ",".join(fi)

    # create DataFrame    
    csites = pd.DataFrame(d)

    # save DataFrame
    csites.to_csv(f"{ROOT}/sites/3cs.txt", sep="\t", index_label="index")
    print(csites)

def get_sstype(seq:str) -> str:

    ''' Get splice site type (based on number and canonicity of NAGs) '''

    # canonical 1-NAG (1C)
    if re.fullmatch("[^N]AG", seq, re.IGNORECASE):
        return "1C"
    
    # non-canonical 1-NAG (1NC)
    if len(seq) == 3:
        return "1NC"

    # canonical NAGNAG (2C)
    if re.fullmatch("[^N]AG[^N]AG", seq, re.IGNORECASE):
        return "2C"
    
    # non-canonical 2-NAG (2C)
    if len(seq) == 6:
        return "2NC"
    
    # 3+ NAG
    if len(seq) % 3 == 0:
        return "3+"

    return "."

# Group cleavage sites
def group_csites():

    ''' Group 3' cleavage sites into 3' splice sites
    
    Source: data/3cs.txt
    
    Destination: data/3ss.txt
    * index (int): 0-indexed row number in DataFrame
    * chrom (str): chromosome on which 3' splice site is found, formatted as 'chr#'
    * strand (str): strand of transcript in which cleavage site is found
        * possible values: "+", "-"
    * num_tri (int): number of cleavage sites + adjacent NAGs in the splice site
    * num_csites (int): number of cleavage sites in the splice site
    * csite_starts (str): comma-separated list of cleavage site down_start positions
    * csite_pos (str): comma-separated list of cleavage site positions (0-indexed) within splice site
    * csite_inds (str): comma-separate list of indices for each cleavage site in 3cs.txt
    * num_scens (str): number of splice scenarios associated with each splice site site
    * ssite_start (int): 0-indexed start position of splice site
        * plus (+) strand: ssite_start is the first base of the splice site motif
        * minus (-) strand: ssite_end is the last base of the upsream intron
    * ssite_end (int): 0-indexed end position of splice site
        * plus (+) strand: ssite_end is the first base of the downstream exon
        * minus (-) strand: ssite_end the last base of the splice site motif
    * ssite_type (str): splice site type, with information about number of trinucleotides and whether each is canonical (NAG)
    * ssite_seq (str): genomic sequence of splice site
    * ssite_iflank_3ss (str): 100 bases upstream of (and not including) splice site
    * ssite_eflank_3ss (str): 100 bases downstream of (and not including) splice site
    * has_nag (int): whether splice site contains a canonical NAG
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
