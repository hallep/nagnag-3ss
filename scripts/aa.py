''' Analyze the effects of NAGNAG alternative splicing on the proteome '''

import pandas as pd
from Bio import Seq
from scripts.sequence import get_aa_transition

# Add amino acid transitions
def add_nagnag_aa_transitions():

    ''' Add variable amino acids to NAGNAG splice scenarios
    
    file: data/nagnag_scens.txt

        "ps_codon" : str
            variable codons when proximally spliced
             * phase -1: "_"
             * phase 0: 3 bases (N2 A G)
             * phase 1: 6 bases (U-1 N2 A G D1 D2)
             * phase 2: 6 bases (U-2 U-1 N A G D1)
        "ds_codon" : str
            variable codons when distally spliced
             * phase -1: "_"
             * phase 0: "."
             * phase 1: 3 bases (U-1 D1 D2)
             * phase 2: 3 bases (U-2 U-1 D1)
        "ps_aa" : str
            variable amino acid(s) when proximally spliced
             * phase -1: "_"
             * phase 0: 1 amino acid
             * phase 1/2: 2 amino acids
        "ds_aa" : str
            variable amino acid when distally spliced
             * phase -1: "_"
             * phase 0: "."
             * phase 1/2: 1 amino acid
    '''

    # load splice scenarios
    scens = pd.read_csv("data/nagnag_scens.txt", sep="\t", index_col=0)

    ps_codon = [""] * len(scens)
    ds_codon = [""] * len(scens)
    ps_aa = [""] * len(scens)
    ds_aa = [""] * len(scens)

    # for each splice scenario:
    for i,(_,scen) in enumerate(scens.iterrows()):

        ps_codon[i], ps_aa[i], ds_codon[i], ds_aa[i] = get_aa_transition(scen["scen_up_seq"], scen["ssite_seq"][3], scen["ssite_down_seq"], scen["frame"])


    # add columns to DataFrame
    scens["ps_codon"] = ps_codon
    scens["ds_codon"] = ds_codon
    scens["ps_aa"] = ps_aa
    scens["ds_aa"] = ds_aa

    # save DataFrame
    scens.to_csv("data/nagnag_scens.txt", sep="\t", index=True)

# Determine possible variable codons
def get_possible_vc():

    ''' Determine all variable codons and their corresponding amino acid transitions

    dst:
     * phase 0: protein/vc_phase0.txt
     * phase 1: protein/vc_phase1.txt
     * phase 2: protein/vc_phase2.txt

        "index" : int
            0-indexed row number in DataFrame
        "ps_codon" : str
            variable codons when proximally spliced
             * phase 0: 3 bases (N2 A G)
             * phase 1: 6 bases (U-1 N2 A G D1 D2)
             * phase 2: 6 bases (U-2 U-1 N A G D1)
        "ds_codon" : str
            variable codons when distally spliced
             * phase 0: "."
             * phase 1: 3 bases (U-1 D1 D2)
             * phase 2: 3 bases (U-2 U-1 D1)
        "ps_aa" : str
            variable amino acid(s) when proximally spliced
             * phase 0: 1 amino acid
             * phase 1/2: 2 amino acids
        "ds_aa" : str
            variable amino acid when distally spliced
             * phase 0: "."
             * phase 1/2: 1 amino acid
    '''

    n_bases = ["A", "C", "G", "T"]

    # PHASE 0
    # get all possible variable codons
    ps_codon_0 = [f"{n}AG" for n in n_bases]

    # get amino acid transitions (translate variable codons)
    ps_aa_0 = [Seq.translate(c) for c in ps_codon_0]

    # create DataFrame
    p0 = pd.DataFrame({
        "ps_codon" : ps_codon_0,
        "ds_codon" : ["."] * 4,
        "ps_aa" : ps_aa_0,
        "ds_aa" : ["."] * 4
    })
    p0 = p0.rename_axis("index").reset_index()
    p0.to_csv("protein/vc_phase0.txt", sep="\t", index=False)

    # PHASE 1
    # get all possible variable codons
    ps_codon_1 = [f"{u1}{n2}AG{d1}{d2}" for u1 in n_bases for n2 in n_bases for d1 in n_bases for d2 in n_bases]
    ds_codon_1 = [f"{u1}{d1}{d2}" for u1 in n_bases for _ in n_bases for d1 in n_bases for d2 in n_bases]

    # get amino acid transitions (translate variable codons)
    ps_aa_1 = [Seq.translate(c) for c in ps_codon_1]
    ds_aa_1 = [Seq.translate(c) for c in ds_codon_1]

    # create & save DataFrame
    p1 = pd.DataFrame({
        "ps_codon" : ps_codon_1,
        "ds_codon" : ds_codon_1,
        "ps_aa" : ps_aa_1,
        "ds_aa" : ds_aa_1
    })
    p1 = p1.rename_axis("index").reset_index()
    p1.to_csv("protein/vc_phase1.txt", sep="\t", index=False)

    # PHASE 2
    # get all possible variable codons
    ps_codon_2 = [f"{u2}{u1}{n2}AG{d1}" for u2 in n_bases for u1 in n_bases for n2 in n_bases for d1 in n_bases]
    ds_codon_2 = [f"{u2}{u1}{d1}" for u2 in n_bases for u1 in n_bases for _ in n_bases for d1 in n_bases]

    # get amino acid transitions (translate variable codons)
    ps_aa_2 = [Seq.translate(c) for c in ps_codon_2]
    ds_aa_2 = [Seq.translate(c) for c in ds_codon_2]

    # create & save DataFrame
    p2 = pd.DataFrame({
        "ps_codon" : ps_codon_2,
        "ds_codon" : ds_codon_2,
        "ps_aa" : ps_aa_2,
        "ds_aa" : ds_aa_2
    })
    p2 = p2.rename_axis("index").reset_index()
    p2.to_csv("protein/vc_phase2.txt", sep="\t", index=False)

# Count variable codons
def count_nagnag_vc(phase:str|int):

    ''' Count the number of NAGNAG splice scenarios associated with each variable codon
    
    src: data/nagnag_scens.txt
    
    dst: protein/vc_phase{phase}.txt

    "num_scen" : int
        number of NAGNAG splice scenarios pertaining to each variable codon
    "num_csp" : int
        number of proximally-spliced splice scenarios in a constitutively-spliced NAGNAG pertaining to each variable codon
    "num_csd" : int
        number of distally-spliced splice scenarios in a constitutively-spliced NAGNAG pertaining to each variable codon
    "num_asp" : int
        number of proximally-spliced splice scenarios in an alternatively-spliced NAGNAG pertaining to each variable codon
    "num_asd" : int
        number of distally-spliced splice scenarios in an alternatively-spliced NAGNAG pertaining to each variable codon
        
    Parameters
    ----------
    phase : str {"0", "1", "2"} or int {0, 1, 2}
        phase/frame of splice scenarios for which to combine variable codons
    '''

    # load splice scenarios
    scens = pd.read_csv("data/nagnag_scens.txt", sep="\t", index_col=0)
    scens = scens[scens["frame"] == int(phase)]

    # load variable codons (by phase)
    vc = pd.read_csv(f"protein/vc_phase{phase}.txt", sep="\t")

    # possible (by motif + cleavage type)
    cp_vc = list(zip(vc["ps_codon"], vc["ds_codon"], ["CSP"] * len(vc)))
    cd_vc = list(zip(vc["ps_codon"], vc["ds_codon"], ["CSD"] * len(vc)))
    ap_vc = list(zip(vc["ps_codon"], vc["ds_codon"], ["ASP"] * len(vc)))
    ad_vc = list(zip(vc["ps_codon"], vc["ds_codon"], ["ASD"] * len(vc)))

    # observed
    obs = list(zip(scens["ps_codon"].str.upper(), scens["ds_codon"].str.upper(), scens["nagnag_ctype"]))

    csp = [obs.count(v) for v in cp_vc]
    csd = [obs.count(v) for v in cd_vc]
    asp = [obs.count(v) for v in ap_vc]
    asd = [obs.count(v) for v in ad_vc]
    num = [cp+cd+ap+ad for cp,cd,ap,ad in zip(csp, csd, asp, asd)]

    # add columns to DataFrame 
    vc["num_scen"] = num
    vc["num_csp"] = csp
    vc["num_csd"] = csd
    vc["num_asp"] = asp
    vc["num_asd"] = asd
    
    # sort
    vc.sort_values(by=["num_scen", "num_asp", "num_csd", "num_asp", "num_asd"],
                   ascending=False, inplace=True)
    
    # save DataFrame
    vc.to_csv(f"protein/vc_phase{phase}.txt", sep="\t", index=False)

# Count amino acid transitions
def count_aa_transitions(phase:str|int):

    ''' Count the number of variable codons per amino acid transition
    
    src: protein/vc_phase{phase}.txt
    
    dst: protein/aat_phase{phase}.txt

        "ps_aa" : str
            variable amino acid(s) when proximally spliced
             * phase 0: 1 amino acid
             * phase 1/2: 2 amino acids
        "ds_aa" : str
            variable amino acid when distally spliced
             * phase 0: "."
             * phase 1/2: 1 amino acid
        "num_vc" : int
            number of variable codons that result in each amino acid transition
        "prop_vc" : int
            proportion variable codons that result in each amino acid transition
        "vc_inds" : str
            comma-separate list of indices for each associated variable codon vc_phase{phase}.txt
        "num_scen" : int
            number of NAGNAG splice scenarios pertaining to each amino acid transition
        "prop_scen" : float
            proportion of NAGNAG splice scenarios (in that phase) pertaining to each amino acid transition
        "num_csp" : int
            number of proximally-spliced splice scenarios in a constitutively-spliced NAGNAG pertaining to each amino acid transition
        "num_csd" : int
            number of distally-spliced splice scenarios in a constitutively-spliced NAGNAG pertaining to each amino acid transition
        "num_asp" : int
            number of proximally-spliced splice scenarios in an alternatively-spliced NAGNAG pertaining to each amino acid transition
        "num_asd" : int
            number of distally-spliced splice scenarios in an alternatively-spliced NAGNAG pertaining to each amino acid transition

    Add aat index to variable codons
        dst: protein/vc_phase{phase}.txt

        "aat_ind" : int
            index of amino acid transition in aat_phase{phase}.txt

    Parameters
    ----------
    phase : str {"0", "1", "2"} or int {0, 1, 2}
        phase/frame of splice scenarios for which to combine variable codons
    '''

    # load variable codons (by phase)
    vc = pd.read_csv(f"protein/vc_phase{phase}.txt", sep="\t", index_col=0)
    vc_unique = list(set(zip(vc["ps_aa"], vc["ds_aa"])))

    aat_inds = [0] * len(vc)

    ps_aa = [c[0] for c in vc_unique]
    ds_aa = [c[1] for c in vc_unique]

    num_vc = [0] * len(vc_unique)
    vc_inds = [""] * len(vc_unique)

    num_scen = [0] * len(vc_unique)
    num_csp = [0] * len(vc_unique)
    num_csd = [0] * len(vc_unique)
    num_asp = [0] * len(vc_unique)
    num_asd = [0] * len(vc_unique)

    # for each unique amino acid transition:
    for i,(ps,ds) in enumerate(zip(ps_aa,ds_aa)):
        these = vc[(vc["ps_aa"] == ps) & (vc["ds_aa"] == ds)]

        num_vc[i] = len(these)
        vc_inds[i] = ",".join(str(ind) for ind in list(these.index))
        num_scen[i] = sum(these["num_scen"])
        num_csp[i] = sum(these["num_csp"])
        num_csd[i] = sum(these["num_csd"])
        num_asp[i] = sum(these["num_asp"])
        num_asd[i] = sum(these["num_asd"])

        # add transition index to variable codons
        for ind in list(these.index):
            aat_inds[ind] = i
    
    aat = pd.DataFrame({
        "ps_aa" : ps_aa,
        "ds_aa" : ds_aa,
        "num_vc" : num_vc,
        "prop_vv" : [n/sum(num_vc) for n in num_vc],
        "vc_inds" : vc_inds,
        "num_scen" : num_scen,
        "prop_scen" : [n/sum(num_scen) for n in num_scen],
        "num_csp" : num_csp,
        "num_csd" : num_csd,
        "num_asp" : num_asp,
        "num_asd" : num_asd
    })
    aat.sort_values(by=["num_scen", "num_asp", "num_csd", "num_asp", "num_asd", "num_vc"],
                    ascending=False, inplace=True)
    aat = aat.rename_axis("index").reset_index()

    # save aat DataFrame
    aat.to_csv(f"protein/aat_phase{phase}.txt", sep="\t", index=False)

    # add "aat_inds" to variable codons
    vc["aat_inds"] = aat_inds
    vc.to_csv(f"protein/vc_phase{phase}.txt", sep="\t", index=True)

# Categorize amino acid transitions
def add_aatype(phase:str|int):

    ''' Categorize each amino acid transition
    
    file: protein/aat_phase{phase}.txt

        "aatype" : str {"DSID", "NSID", "CSID", "IDR"}
            amino acid transition type
        "DSID": duplicate simple indel
            p1 == d; p2 == d
        "NSID": N-terminus (pre) simple indel
            p1 != d; p2 == d
        "CSID": C-terminus (post) simple indel
            p1 == d; p2 != d
        "IDR": indel + replacement
            p1 != d; p2 != d

    Parameters
    ----------
    phase : str {"1", "2"} or int {1, 2}
        phase/frame of splice scenarios for which to combine variable codons
    '''

    # load amino acid transitions (by phase)
    aat = pd.read_csv(f"protein/aat_phase{phase}.txt", sep="\t", index_col=0)

    aatype = [""] * len(aat)

    # for each transition:
    for i,(_,t) in enumerate(aat.iterrows()):
        p1 = t["ps_aa"][0]
        p2 = t["ps_aa"][1]
        d = t["ds_aa"]

        # duplicate simple indel
        if (p1 == d) and (p2 == d):
            aatype[i] = "DSID"

        # C-terminus (post) simple indel
        elif (p1 == d):
            aatype[i] = "CSID"

        # N-terminus (pre) simple indel
        elif (p2 == d):
            aatype[i] = "NSID"

        # indel + replacement
        else:
            aatype[i] = "IDR"

    # add column to DataFrame
    aat["aatype"] = aatype

    # save DataFrame
    aat.to_csv(f"protein/aat_phase{phase}.txt", sep="\t", index=True)

# Count amino acid transition types
def count_aatype(phase:str|int):

    ''' Count the number of each amino acid transition type
    
    src: protein/aat_phase{phase}.txt

    dst: protein/aatype_phase{phase}.txt

        "aatype" : str {"DSID", "NSID", "CSID", "IDR"}
            amino acid transition type
             * "DSID": duplicate simple indel (p1 == d; p2 == d)
             * "NSID": N-terminus (pre) simple indel (p1 != d; p2 == d)
             * "CSID": C-terminus (post) simple indel (p1 == d; p2 != d)
             * "IDR": indel + replacement (p1 != d; p2 != d)
        "num_vc" : int
            number of variable codons that result in each amino acid transition type
        "prop_vc" : int
            proportion variable codons that result in each amino acid transition type
        "num_aat" : int
            number of amino acid transitions in each type
        "prop_aat" : int
            proportion of amino acid transitions in each type
        "num_scen" : int
            number of NAGNAG splice scenarios pertaining to each amino acid transition type
        "prop_scen" : float
            proportion of NAGNAG splice scenarios (in that phase) pertaining to each amino acid transition type
        "oe_ratio" : float
            ratio of observed / expected
             * observed: prop_scen
             * expected: prop_vc

    Parameters
    ----------
    phase : str {"1", "2"} or int {1, 2}
        phase/frame of splice scenarios for which to combine variable codons
    '''

    # load amino acid transitions
    aat = pd.read_csv(f"protein/aat_phase{phase}.txt", sep="\t", index_col=0)

    aatype = ["DSID", "NSID", "CSID", "IDR"]

    num_vc = [0] * len(aatype)
    num_aat = [0] * len(aatype)
    num_scen = [0] * len(aatype)

    # for each amino acid transition type
    for i,aaty in enumerate(aatype):
        these = aat[aat["aatype"] == aaty]

        num_vc[i] = sum(these["num_vc"])
        num_aat[i] = len(these)
        num_scen[i] = sum(these["num_scen"])
    
    # get proportions
    prop_vc = [n/sum(num_vc) for n in num_vc]
    prop_aat = [n/sum(num_aat) for n in num_aat]
    prop_scen = [n/sum(num_scen) for n in num_scen]

    # get ratios
    oe_ratio = [o/e for o,e in zip(prop_scen, prop_vc)]

    # create DataFrame
    aatypes = pd.DataFrame({
        "aatype" : aatype,
        "num_vc" : num_vc,
        "prop_vc" : prop_vc,
        "num_aat" : num_aat,
        "prop_aat" : prop_aat,
        "num_scen" : num_scen,
        "prop_scen" : prop_scen,
        "oe_ratio" : oe_ratio
    })

    # save DataFrame
    aatypes.to_csv(f"protein/aatype_phase{phase}.txt", sep="\t", index=False)
