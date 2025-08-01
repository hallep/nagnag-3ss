''' Get frequencies:
* splice site types
* 1-NAG motifs
* NAGNAG splice types
* NAGNAG splice scenario phases
'''

import pandas as pd

# Splice site types
def count_ss_types():

    ''' Count the number of each uppercase splice site type
    
    src: data/3ss.txt
    
    dst: stats/sstype_count.txt

        "sstype": str {"1C", "1NC", "2C", "2NC", "3+"}
            splice site type
        "num_ss": int
            number of splice sites of each type
        "prop_ss": float
            proportion of splice sites of each type
        "num_cs": int
            number of cleavage sites of each type
        "prop_ss": float
            proportion of cleavage sites of each type
        "num_scen": int
            number of splice scenarios of each type
        "prop_scen" : float
            proportion of splice scenarios of each type
    '''
    
    # load splice sites
    ssites = pd.read_csv("data/3ss.txt", sep="\t", index_col=0)
    ssites = ssites[ssites["uppercase"] == 1]

    sstypes = ["1C", "1NC", "2C", "2NC", "3+"]

    # get counts
    count_ss = [len(ssites[ssites["ssite_type"] == s]) for s in sstypes]
    count_cs = [sum(ssites[ssites["ssite_type"] == s]["num_csites"]) for s in sstypes]
    count_scen = [sum(ssites[ssites["ssite_type"] == s]["num_scens"]) for s in sstypes]

    # create DataFrame
    df = pd.DataFrame({
        "sstype" : sstypes + ["total"],
        "num_ss" : count_ss + [sum(count_ss)],
        "prop_ss" : [c/sum(count_ss) for c in count_ss] + [1],
        "num_cs" : count_cs + [sum(count_cs)],
        "prop_cs" : [c/sum(count_cs) for c in count_cs] + [1],
        "num_scen" : count_scen + [sum(count_scen)],
        "prop_scen" : [c/sum(count_scen) for c in count_scen] + [1]
    })

    # save DataFrame
    df.to_csv("stats/sstype_count.txt", sep="\t", index=False)

    print("=== 3' Splice Sites ===")
    print(df)

# 1-NAG motifs
def count_1nag_motifs():

    ''' Count the number of uppercase canonical NAGs of each motif
    
    src: data/3ss.txt
    
    dst: stats/canon_nag_motif_count.txt

        "motif" : str {"AAG", "CAG", "GAG", "TAG"}
            1-NAG motif
        "num_ss" : int
            number of 1-NAG splice sites of each motif
        "prop_ss" : float
            proportion of 1-NAG splice sites of each motif
        "num_cs" : int
            number of 1-NAG cleavage sites of each motif
        "prop_ss" : float
            proportion of 1-NAG cleavage sites of each motif
        "num_scen" : int
            number of 1-NAG splice scenarios of each motif
        "prop_scen" : float
            proportion of 1-NAG splice scenarios of each motif
    '''

    # load splice sites
    ssites = pd.read_csv("data/3ss.txt", sep="\t", index_col=0)
    ssites = ssites[ssites["uppercase"] == 1]

    motifs = ["AAG", "CAG", "GAG", "TAG"]

    # get counts
    count_ss = [len(ssites[(ssites["ssite_type"] == "1C") & (ssites["ssite_seq"] == m)]) for m in motifs]
    count_cs = [sum(ssites[(ssites["ssite_type"] == "1C") & (ssites["ssite_seq"] == m)]["num_csites"]) for m in motifs]
    count_scen = [sum(ssites[(ssites["ssite_type"] == "1C") & (ssites["ssite_seq"] == m)]["num_scens"]) for m in motifs]

    # create DataFrame
    df = pd.DataFrame({
        "motif" : motifs + ["total"],
        "num_ss" : count_ss + [sum(count_ss)],
        "prop_ss" : [c/sum(count_ss) for c in count_ss] + [1],
        "num_cs" : count_cs + [sum(count_cs)],
        "prop_cs" : [c/sum(count_cs) for c in count_cs] + [1],
        "num_scen" : count_scen + [sum(count_scen)],
        "prop_scen" : [c/sum(count_scen) for c in count_scen] + [1]
    })

    # save DataFrame
    df.to_csv("stats/canon_nag_motif_count.txt", sep="\t", index=False)

    print("\n=== 1-NAG Motifs ===")
    print(df)

# NAGNAG splice types (PS, DS, AS)
def count_nagnag_splice_types():
    
    '''Count the number of uppercase canonical NAGNAGs of each splice type
    
    src: data/3ss.txt
    
    dst: stats/nagnag_stype_count.txt

        "stype" : str {"PS", "DS", "AS"}
            NAGNAG splice type
        "num_ss" : int
            number of NAGNAG splice sites of each splice type
        "prop_ss" : float
            proportion of NAGNAG splice sites of each splice type
        "num_cs" : int
            number of NAGNAG cleavage sites of each splice type
        "prop_ss" : float
            proportion of NAGNAG cleavage sites of each splice type
        "num_scen" : int
            number of NAGNAG splice scenarios of each splice type
        "prop_scen" : float
            proportion of NAGNAG splice scenarios of each splice type
    '''

    # load splice sites
    ssites = pd.read_csv("data/3ss.txt", sep="\t", index_col=0)
    ssites = ssites[ssites["uppercase"] == 1]

    stype = ["PS", "DS", "AS"]
    pos = ["0", "1", "0,1"]

    # get counts
    count_ss = [len(ssites[(ssites["ssite_type"] == "2C") & (ssites["csite_pos"] == p)]) for p in pos]
    count_cs = [sum(ssites[(ssites["ssite_type"] == "2C") & (ssites["csite_pos"] == p)]["num_csites"]) for p in pos]
    count_scen = [sum(ssites[(ssites["ssite_type"] == "2C") & (ssites["csite_pos"] == p)]["num_scens"]) for p in pos]

    # create DataFrame
    df = pd.DataFrame({
        "stype" : stype + ["total"],
        "num_ss" : count_ss + [sum(count_ss)],
        "prop_ss" : [c/sum(count_ss) for c in count_ss] + [1],
        "num_cs" : count_cs + [sum(count_cs)],
        "prop_cs" : [c/sum(count_cs) for c in count_cs] + [1],
        "num_scen" : count_scen + [sum(count_scen)],
        "prop_scen" : [c/sum(count_scen) for c in count_scen] + [1]
    })

    # save DataFrame
    df.to_csv("stats/nagnag_stype_count.txt", sep="\t", index=False)

    print(f"\n=== NAGNAG Splice Types ===")
    print(df)

# NAGNAG splice scenario phases
def count_nagnag_phases():
    
    '''Count the number of NAGNAG splice scenarios of each phase
    
    src: data/nagnag_scens.txt
    
    dst: stats/nagnag_scen_phase_count.txt

        "phase" : int {-1, 0, 1, 2}
            phase of downstream exon
             * -1: non-CDS
        "count" : int
            number of NAGNAG splice scenarios in each phase
        "prop" : float
            proportion of NAGNAG splice scenarios in each phase
        "prop_CDS" : float
            proportion of NAGNAG splice scenarios in each CDS phase
    '''

    # load splice scenarios
    scens = pd.read_csv("data/nagnag_scens.txt", sep="\t", index_col=0)

    frame = [-1, 0, 1, 2]
    
    # get counts
    count = [len(scens[scens["frame"] == f]) for f in frame]

    # create DataFrame
    df = pd.DataFrame({
        "phase" : frame + ["total"],
        "count" : count + [sum(count)],
        "prop" : [c/sum(count) for c in count] + [1],
        "prop_CDS" : [0] + [c/sum(count[1:]) for c in count[1:]] + [1]
    })

    # save DataFrame
    df.to_csv("stats/nagnag_scen_phase_count.txt", sep="\t", index=False)

    print(f"\n=== NAGNAG Splice Scenario Phases ===")
    print(df)
