from tqdm import tqdm # type: ignore
import pandas as pd # type: ignore
from sequence import get_aa_transition

''' NAGNAG SPLICE SCENARIOS '''

# GET NAGNAG SPLICE INSTANCES
def get_splice_instances():
    
    ''' Get NAGNAG splice scenarios (from thesis) as splice instances
    
    src: /mnt/data_3/hallep/nagnag/thesis/nagnag_scens.txt

    dst: /mnt/data_3/hallep/nagnag/nagnag_splice_instances.txt
    '''

    insts = pd.read_csv("/mnt/data_3/hallep/nagnag/thesis/nagnag_scens.txt", sep="\t", index_col=0)
    insts.to_csv("/mnt/data_3/hallep/nagnag/nagnag_splice_instances.txt", sep="\t", index=True)

# GROUP NAGNAG SPLICE INSTANCES
def group_splice_instances():

    ''' Group NAGNAG splice instances

    src:
    * splice sites: /mnt/data_3/hallep/nagnag/nagnag_3ss.txt
    * splice instances: /mnt/data_3/hallep/nagnag/nagnag_splice_instances.txt
    
    dst:
    * all: /mnt/data_3/hallep/nagnag/nagnag_splice_scenarios.txt
    * cds: /mnt/data_3/hallep/nagnag/nagnag_splice_scenarios_cds.txt

    Columns:
    * "index" : int \\
        0-indexed row number in (all) DataFrame
    * "accession" : str \\
        NCBI accession numbers, split by "/" for the same cleavage site and "," for different cleavage sites
    * "chromosome" : str \\
        chromosome, in the format "chr#"
    * "downstream_start" : int \\
        0-indexed position
    * "strand" : str {"+", "-"} \\
        strand of transcript
    * "upstream_end" : int \\
        0-indexed position of first intronic base after upstream exon
    * "scen_splice_type" : str {"AS", "PS", "DS"} \\
        splice type of splice scenario
        * "AS" : there are two splice scenarios corresponding to the same upstream exon, phase, and rtype
        * "PS" : there is one splice scenario for the proximal NAG
        * "DS" : there is one splice scenario for the distal NAG
    * "nagnag_splice_type" : str {"AS", "PS", "DS"} \\
        splice type of NAGNAG
        * "AS" : both cleavage sites of NAGNAG are used
        * "PS" : only the first NAG is annotated as being used
        * "DS" : only the distal NAG is annotated as being used
    * "rtype" : str {"CDS", "ncRNA", "5'UTR", "3'UTR"} \\
        RNA type
        * "CDS" : in coding sequence of a coding RNA
        * "ncRNA" : in a non-coding RNA
        * "5'UTR" : in the 5' UTR of a coding RNA
        * "3'UTR" : in the 3' UTR of a coding RNA
    * "phase" : int {-1, 0, 1, 2} \\
        phase of downstream exon
    * "u-2" : str \\
        second-to-last base of upstream exon
    * "u-1" : str \\
        last base of upstream exon
    * "n1" : str \\
        identity of first NAG
    * "n2" : str \\
        identity of second NAG
    * "d1" : str \\
        first base after NAGNAG
    * "d2" : str \\
        second base after NAGNAG
    * "ps_vc" : str \\
        variable codons when proximally spliced
    * "ds_vc" : str \\
        variable codons when distally spliced
    * "ps" : str \\
        amino acids when proximally spliced
    * "ds" : str \\
        amino acids when distally spliced
    * "ps_A" : str \\
        amino acids when proximally spliced + n2 = A
    * "ps_C" : str \\
        amino acids when proximally spliced + n2 = C
    * "ps_G" : str \\
        amino acids when proximally spliced + n2 = G
    * "ps_T" : str \\
        amino acids when proximally spliced + n2 = T
    * "ssite_ind" : int \\
        index in /mnt/data_3/hallep/nagnag/nagnag_3ss.txt of associated NAGNAG splice site
    * "sscen_ind" : str \\
        index in /mnt/data_3/hallep/nagnag/nagnag_splice_instances.txt of associated NAGNAG splice scenario(s) \\
        if alternatively spliced: {index of proximal},{index of distal}
    '''

    ssites = pd.read_csv("/mnt/data_3/hallep/nagnag/nagnag_3ss.txt", sep="\t", index_col=0)
    sinsts = pd.read_csv("/mnt/data_3/hallep/nagnag/nagnag_splice_instances.txt", sep="\t", index_col=0)

    access = []
    chrom = []
    down_start = []
    strand = []
    up_end = []
    scen_stype = []
    site_stype = []
    rtype = []
    phase = []
    u2 = []
    u1 = []
    n1 = []
    n2 = []
    d1 = []
    d2 = []
    ssite_ind = []
    sscen_ind = []

    # for each splice site:
    for s,ssite in ssites.iterrows():

        # get splice instances
        inds = [int(j) for i in ssite["scen_inds"].split(",") for j in i.split("/")]
        insts = sinsts.loc[inds]

        # important variable components: upstream end, frame, RNA type, processed {0,1}
        scen_ids = [list(z) for z in zip(insts["up_end"], insts["scen_up_seq"], insts["ssite_down_seq"], insts["frame"], insts["rtype"], [0] * len(insts))]

        # for each splice instances:
        for i,(x,inst) in enumerate(insts.iterrows()):

            # if unprocessed:
            if scen_ids[i][5] == 0:

                # mark as processed            
                scen_ids[i][5] = 1

                useq = inst["scen_up_seq"]
                sseq = inst["ssite_seq"]
                dseq = inst["ssite_down_seq"]

                # search for NAGNAG alternative splicing
                id = [inst["up_end"], inst["scen_up_seq"], inst["ssite_down_seq"], inst["frame"], inst["rtype"], 0]

                # alternatively spliced
                if id in scen_ids:
                    
                    j = scen_ids.index(id)
                    y = list(insts.index)[j]

                    # mark as processed
                    scen_ids[j][5] = 1

                    if inst["nagnag_ctype"][-1] == "P" and insts["nagnag_ctype"][y][-1] == "D":
                        access.append(f"{inst["accession"]},{sinsts["accession"][y]}")
                        sscen_ind.append(f"{x},{y}")
                    elif inst["nagnag_ctype"][-1] == "D" and insts["nagnag_ctype"][y][-1] == "P":
                        access.append(f"{sinsts["accession"][y]},{inst["accession"]}")
                        sscen_ind.append(f"{y},{x}")

                    scen_stype.append("AS")

                # constitutively spliced
                else:
                    
                    access.append(inst["accession"])
                    sscen_ind.append(str(x))

                    # proximally spliced
                    if inst["nagnag_ctype"][-1] == "P":
                        scen_stype.append("PS")

                    # distally spliced
                    elif inst["nagnag_ctype"][-1] == "D":
                        scen_stype.append("DS")
                                    
                # NAGNAG type
                if ssite["csite_pos"] == "0,1": 
                    site_stype.append("AS")
                elif ssite["csite_pos"] == "0":
                    site_stype.append("PS")
                elif ssite["csite_pos"] == "1":
                    site_stype.append("DS")

                chrom.append(inst["chrom"])
                down_start.append(inst["ssite_end"])
                strand.append(inst["strand"])
                up_end.append(inst["up_end"])
                rtype.append(inst["rtype"])
                phase.append(inst["frame"])
                u2.append(useq[-2])
                u1.append(useq[-1])
                n1.append(sseq[0])
                n2.append(sseq[3])
                d1.append(dseq[0])
                d2.append(dseq[1])

                ssite_ind.append(s)

    # create DataFrame
    df = pd.DataFrame({
        "accession" : access,
        "chromosome" : chrom,
        "downstream_start" : down_start,
        "strand" : strand,
        "upstream_end" : up_end,
        "scen_splice_type" : scen_stype,
        "nagnag_splice_type" : site_stype,
        "rtype" : rtype,
        "phase" : phase,
        "u-2" : u2,
        "u-1" : u1,
        "n1" : n1,
        "n2" : n2,
        "d1" : d1,
        "d2" : d2,
        "ps_vc" : [get_aa_transition(f"N{u_2}{u_1}", n_2, f"{d_1}{d_2}N", p)[0] for u_2,u_1,n_2,d_1,d_2,p in zip(u2,u1,n2,d1,d2,phase)],
        "ds_vc" : [get_aa_transition(f"N{u_2}{u_1}", n_2, f"{d_1}{d_2}N", p)[1] for u_2,u_1,n_2,d_1,d_2,p in zip(u2,u1,n2,d1,d2,phase)],
        "ps" : [get_aa_transition(f"N{u_2}{u_1}", n_2, f"{d_1}{d_2}N", p)[2] for u_2,u_1,n_2,d_1,d_2,p in zip(u2,u1,n2,d1,d2,phase)],
        "ds" : [get_aa_transition(f"N{u_2}{u_1}", n_2, f"{d_1}{d_2}N", p)[3] for u_2,u_1,n_2,d_1,d_2,p in zip(u2,u1,n2,d1,d2,phase)],
        "ps_A" : [get_aa_transition(f"N{u_2}{u_1}", "A", f"{d_1}{d_2}N", p)[2] for u_2,u_1,d_1,d_2,p in zip(u2,u1,d1,d2,phase)],
        "ps_C" : [get_aa_transition(f"N{u_2}{u_1}", "C", f"{d_1}{d_2}N", p)[2] for u_2,u_1,d_1,d_2,p in zip(u2,u1,d1,d2,phase)],
        "ps_G" : [get_aa_transition(f"N{u_2}{u_1}", "G", f"{d_1}{d_2}N", p)[2] for u_2,u_1,d_1,d_2,p in zip(u2,u1,d1,d2,phase)],
        "ps_T" : [get_aa_transition(f"N{u_2}{u_1}", "T", f"{d_1}{d_2}N", p)[2] for u_2,u_1,d_1,d_2,p in zip(u2,u1,d1,d2,phase)],
        "ssite_ind" : ssite_ind,
        "sscen_ind" : sscen_ind
    })
    df = df.rename_axis("index").reset_index()

    # save all scenarios
    df.to_csv("/mnt/data_3/hallep/nagnag/nagnag_splice_scenarios.txt", sep="\t", index=False)

    # save CDS scenarios
    df[df["rtype"] == "CDS"].to_csv("/mnt/data_3/hallep/nagnag/nagnag_splice_scenarios_cds.txt", sep="\t", index=False)

    print(df)

# get_splice_instances()
# group_splice_instances()

''' 1-NAG SPLICE SCENARIOS '''

# GET 1-NAG CLEAVAGE SITES
def get_1nag_csites():

    ''' Isolate 1-NAG cleavage sites

    src:
    * splice sites: /mnt/data_3/hallep/nagnag/3ss.txt \\
    * cleavage sites: /mnt/data_3/hallep/nagnag/3cs.txt
    
    dst: /mnt/data_3/hallep/nagnag/1nag_3cs.txt
    '''
    
    # get 1-NAG splice sites
    ssites = pd.read_csv("/mnt/data_3/hallep/nagnag/3ss.txt", sep="\t", index_col=0)
    ssites = ssites[(ssites["uppercase"] == 1) & (ssites["ssite_type"] == "1C")]

    # get 1-NAG cleavage sites
    csites = pd.read_csv("/mnt/data_3/hallep/nagnag/3cs.txt", sep="\t", index_col=0)
    cs_inds = [int(i) for i in ssites["csite_inds"]]
    csites = csites.loc[cs_inds].sort_index()

    # save 1-NAG cleavage sites
    csites.to_csv("/mnt/data_3/hallep/nagnag/1nag_3cs.txt", sep="\t", index=True)
    print(csites)

# SPLIT 1-NAG SPLICE SCENARIOS
def split_1nag_scens():

    ''' Split 1-NAG cleavage sites into scenarios

    cleavage sites: /mnt/data_3/hallep/nagnag/1nag_3cs.txt
    * "scen_inds" : str \\
        comma-separated list of indices in 1nag_splice_scenarios.txt
    
    splice scenarios: /mnt/data_3/hallep/nagnag/1nag_splice_scenarios.txt
    * "index" : int \\
        1-index of 
    * "chromosome" : str \\
        chromosome, in the format "chr#"
    * "downstream_start" : int \\
        0-indexed position
    * "strand" : str {"+", "-"} \\
        strand of transcript
    * "rtype" : str {"CDS", "ncRNA", "5'UTR", "3'UTR"} \\
        RNA type
        * "CDS" : in coding sequence of a coding RNA
        * "ncRNA" : in a non-coding RNA
        * "5'UTR" : in the 5' UTR of a coding RNA
        * "3'UTR" : in the 3' UTR of a coding RNA
    * "upstream_end" : int \\
        0-indexed position of first intronic base after upstream exon
    * "phase" : int {-1, 0, 1, 2} \\
        phase of downstream exon
    * "accession" : str \\
        NCBI accession numbers, split by "/" for the same cleavage site and "," for different cleavage sites
    * "csite_seq" : str {"CAG", "TAG", "AAG", "GAG"} \\
        cleavage site sequence
    * "eflank_5ss" : str \\
        last 100 bases of the upstream exon
    * "iflank_3ss" : str \\
        first 100 bases of the intron following the upstream exon
    * "iflank_3ss" : str \\
        last 100 bases of the intron (including csite_seq)
    * "eflank_3ss" : str \\
        first 100 bases of the downstream exon

    * "csite_ind" : int \\
        index in 3cs.txt of associated splice cleavage site
        
    '''

    # load cleavage sites
    csites = pd.read_csv("/mnt/data_3/hallep/nagnag/1nag_3cs.txt", sep="\t", index_col=0)
    scen_ind = [[] for _ in range(len(csites))]

    chrom = []
    down_start = []
    strand = []
    rtype = []
    up_end = []
    phase = []
    accession = []
    csite_seq = []
    eflank_5ss = []
    iflank_5ss = []
    iflank_3ss = []
    eflank_3ss = []
    cs_ind = []

    x = 0
    # for each cleavage site:
    for c,(_,cs) in tqdm(enumerate(csites.iterrows())):
        
        rtI = cs["rtype"].split(",")
        upI = cs["up_end"].split(",")
        phI = cs["frame"].split(",")
        acI = cs["accession"].split(",")
        efI = cs["scen_eflank_5ss"].split(",")
        ifI = cs["scen_iflank_5ss"].split(",")
        
        # for each splice scenario:
        for r,u,p,a,e,i in zip(rtI, upI, phI, acI, efI, ifI):
            
            # scenario details
            chrom.append(cs["chrom"])
            down_start.append(cs["down_start"])
            strand.append(cs["strand"])
            rtype.append(r)
            up_end.append(u)
            phase.append(p)
            accession.append(a)
            csite_seq.append(cs["csite_seq"])
            eflank_5ss.append(e)
            iflank_5ss.append(i)
            iflank_3ss.append(cs["csite_iflank_3ss"])
            eflank_3ss.append(cs["csite_eflank_3ss"])
            cs_ind.append(c)
            
            # scenario index
            scen_ind[c].append(str(x))
            x += 1

    # create DataFrame    
    df = pd.DataFrame({
        "chrom" : chrom,
        "down_start" : down_start,
        "strand" : strand,
        "rtype" : rtype,
        "up_end" : up_end,
        "phase" : phase,
        "accession" : accession,
        "csite_seq" : csite_seq,
        "eflank_5ss" : eflank_5ss,
        "iflank_5ss" : iflank_5ss,
        "iflank_3ss" : iflank_3ss,
        "eflank_3ss" : eflank_3ss,
        "csite_ind" : cs_ind
    })

    # save splice scenarios
    df.to_csv("/mnt/data_3/hallep/nagnag/1nag_splice_scenarios.txt", sep="\t", index=True, index_label="index")

    print("\n===== SPLICE SCENARIOS =====")
    print(df)

    # save cleavage sites
    csites["scen_inds"] = [",".join(i) for i in scen_ind]
    csites.to_csv("/mnt/data_3/hallep/nagnag/1nag_3cs.txt", sep="\t", index=True)

    print("\n===== CLEAVAGE SITES =====")
    print(csites)

# get_1nag_csites()
# split_1nag_scens()

''' .BED FILE '''

# CREATE NAGNAG .BED FILE
def nagnag_bed():

    ''' Create .bed file for NAGNAG splice scenarios
    
    src: /mnt/data_3/hallep/nagnag/nagnag_splice_scenarios.txt

    dst: /mnt/data_3/hallep/reference/tx_data/nagnags.bed
    * name: {index}_{motif}_{rtype}_phase{phase}_{strand}_{stype}
    '''

    colors = {
        "AS" : {
            "CDS" : "4,97,152",
            "5'UTR" : "117,159,194",
            "3'UTR" : "117,159,194",
            "ncRNA" : "1,45,74"
        },
        "PS" : {
            "CDS" : "0,96,52",
            "5'UTR" : "114,157,128",
            "3'UTR" : "114,157,128",
            "ncRNA" : "0,44,21",
        },
        "DS" : {
            "CDS" : "132,0,122",
            "5'UTR" : "185,117,175",
            "3'UTR" : "185,117,175",
            "ncRNA" : "63,0,58"
        }
    }

    # load splice scenarios
    scens = pd.read_csv("/mnt/data_3/hallep/nagnag/nagnag_splice_scenarios.txt", sep="\t", index_col=0)

    start = [0] * len(scens)
    end = [0] * len(scens)
    name = [""] * len(scens)
    tStart = [0] * len(scens)
    tEnd = [0] * len(scens)
    rgb = [""] * len(scens)

    # for each splice scenario:
    for i,(x,scen) in enumerate(scens.iterrows()):
        
        if scen["strand"] == "+":
            start[i] = scen["downstream_start"] - 6
            end[i] = scen["downstream_start"]
            tStart[i] = scen["downstream_start"] - 6 if scen["rtype"] == "CDS" else "."
            tEnd[i] = scen["downstream_start"] if scen["rtype"] == "CDS" else "."
        
        elif scen["strand"] == "-":
            start[i] = scen["downstream_start"]
            end[i] = scen["downstream_start"] + 6
            tStart[i] = scen["downstream_start"] if scen["rtype"] == "CDS" else "."
            tEnd[i] = scen["downstream_start"] + 6 if scen["rtype"] == "CDS" else "."

        name[i] = f"{x}_{scen["n1"]}AG{scen["n2"]}AG_{scen["rtype"]}_phase{scen["phase"]}_{scen["strand"]}_{scen["scen_splice_type"]}"

        rgb[i] = colors[scen["scen_splice_type"]][scen["rtype"]]

    # create .bed file
    bed = pd.DataFrame({
        "chrom" : scens["chromosome"].values,
        "start" : start,
        "end" : end,
        "name" : name,
        "score" : [0] * len(scens),
        "strand" : scens["strand"].values,
        "thickStart" : tStart,
        "thickEnd" : tEnd,
        "rgb" : rgb,
    })
    bed.to_csv("/mnt/data_3/hallep/reference/tx_data/nagnags.bed", sep="\t", index=False, header=False)

    print(bed)

# nagnag_bed()
