import itertools
from tqdm import tqdm
import pandas as pd
import numpy as np
import scipy.stats as stats
from sequence import get_transition, get_binary_outcomes, get_categorical_outcomes, TSNS, N, AA, BIN, CAT
from snps_tools import CHROMS, dtypes, var_info, snp_dbs, coord_col, maf_col, nn_effects

dtypes.update({
    "alt_ind" : "int",
    "ref_single" : "int",
    "alt_single" : "int",
})

''' ===== SNP Events & Scenarios ===== '''

# SNP Events
def get_snp_events(db:str):
    
    ''' Get SNP events

    source:
    * cleavage sites: /mnt/data_3/hallep/nagnag/3cs.txt
    * sites: /mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_containing_ssites.txt
    * SNPs: /mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snps.txt

    destination:
    * /mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_events.txt
        * columns: "index", "chrom", "motif_start", "motif_end", "strand", "splice_type",
            "ref_seq", "alt_seq", "down_seq", "vnt_id", "vnt_coord", "vnt_pos", "ref", "alt", "vnt_effect",
            "vsite_ind", "nsite_ind", "ssite_ind", "csite_pos", "csite_inds", {maf_col[db][dst]}
        * "index" is new
    * /mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snps.txt
        * added column: "event_inds"

    Parameters
    ----------
    db : str {"dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice}
        SNP database
    '''

    # load SNP-containing sites
    sites = pd.read_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_containing_ssites.txt", sep="\t", index_col=0, dtype=dtypes)

    # load NAGNAG-affecting SNPs
    snps = pd.read_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snps.txt", sep="\t", index_col=0, dtype=dtypes)

    # SNP scenario data
    d = {
        "chrom" : [],
        "motif_start" : [],
        "motif_end" : [],
        "strand" : [],
        "splice_type" : [],
        "ref_seq" : [],
        "alt_seq" : [],
        "down_seq" : [],
        "vnt_id" : [],
        "vnt_coord" : [],
        "vnt_pos" : [],
        "ref" : [],
        "alt" : [],
        "vnt_effect" : [],
    }
    d.update({c : [] for c in maf_col[db]["dst"]})
    d.update({
        "vsite_ind" : [],
        "nsite_ind" : [],
        "ssite_ind" : [],
        "csite_pos" : [],
        "csite_inds" : [],
    })
    
    # for each SNP:
    for i,snp in snps.iterrows():

        # alternative allele
        alt = snp["alt"].split(",")[int(snp["alt_ind"])]

        v = len(snp["vsite_ind"].split(","))
        n = len(snp["nsite_ind"].split(","))
        s = len(snp["ssite_ind"].split(","))

        if (v != n) or (v != s) or (n != s):
            print(snp["vsite_ind"], snp["nsite_ind"], snp["ssite_ind"])

        # site
        for v,n,s in zip(snp["vsite_ind"].split(","), snp["nsite_ind"].split(","), snp["ssite_ind"].split(",")):
            site = sites.loc[v]
            st = "PS" if site["csite_pos"] == "0" else "DS" if site["csite_pos"] == "1" else "AS"
            
            # alternative allele (by strand)
            rc = {"A":"T", "T":"A", "C":"G", "G":"C"}
            r = snp["ref"] if (var_info[db][4]) or (site["strand"] == "+") else rc[snp["ref"]]
            a = alt if (site["strand"] == "+") or (db == "HGMD_splice") else rc[alt]

            # SNP sequence
            pos = {n : int(p) for n,p in zip(site["snp_ids"].split(","), site["snp_pos"].split(","))}[i]
            motif = site["motif_seq"][:pos] + a + site["motif_seq"][pos+1:]
            
            # location information (first base after motif)
            d["chrom"].append(site["chrom"])
            d["motif_start"].append(site["motif_start"])
            d["motif_end"].append(site["motif_end"])
            d["strand"].append(site["strand"])
            d["splice_type"].append(st)
            
            # sequences
            d["ref_seq"].append(site["motif_seq"])
            d["alt_seq"].append(motif)
            d["down_seq"].append(site["motif_eflank_3ss"][:3])

            # variant info: ID, coordinate, position, ref, alt, type
            d["vnt_id"].append(i)
            d["vnt_coord"].append(snp[coord_col[db]] - 1)
            d["vnt_pos"].append(pos)
            d["ref"].append(r)
            d["alt"].append(a)
            d["vnt_effect"].append(snp["snp_type"])

            # reference indices
            d["vsite_ind"].append(v)
            d["nsite_ind"].append(n)
            d["ssite_ind"].append(s)
            d["csite_pos"].append(site["csite_pos"])
            d["csite_inds"].append(site["csite_inds"])

            # minor allele frequencies
            for src,dst in zip(maf_col[db]["src"], maf_col[db]["dst"]):
                d[dst].append(snp[src])

    # DataFrame
    events = pd.DataFrame(d)

    # sort
    events["order"] = events["chrom"].apply(lambda x: CHROMS.index(x[3:]))
    events.sort_values(by=["order", "motif_start", "strand", "vnt_pos"], ignore_index=True, inplace=True)
    events.drop(columns="order", inplace=True)

    events.to_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_events.txt", sep="\t", index_label="index")

    print("\n===== SNP Events =====")
    print(events)

    # add scenario indices to SNPs
    snps["event_inds"] = [",".join(map(str, events[events["vnt_id"] == i].index.values)) for i in snps.index]
    snps.to_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snps.txt", sep="\t")

# SNP Scenarios
def get_snp_scenarios(db:str):
    
    ''' Get SNP scenarios

    source:
    * cleavage sites: /mnt/data_3/hallep/nagnag/3cs.txt
    * SNP events: /mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_events.txt

    destination:
    * /mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_scenarios.txt
        * columns: "index", "chrom", "down_start", "strand", "up_end", "rtype", "phase", "splice_type",
            "up_seq", "ref_seq", "alt_seq", "down_seq",
            "vnt_id", "vnt_coord", "vnt_pos", "ref", "alt", "vnt_effect",
            {maf columns}, "vsite_ind", "ssite_ind", "csite_ind", "event_ind"
        * "index" is new
    * /mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snps.txt
        * added column: "scen_inds"
    * /mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_events.txt
        * added column: "scen_inds"

    Parameters
    ----------
    db : str {"dbSNP", "ClinVar", "HGMD_splice}
        SNP database
    '''

    # load cleavage sites
    csites = pd.read_csv("/mnt/data_3/hallep/nagnag/3cs.txt", sep="\t", index_col=0, dtype=dtypes)

    # load NAGNAG-affecting SNPs
    snps = pd.read_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snps.txt", sep="\t", index_col=0, dtype=dtypes)

    # load SNP events
    events = pd.read_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_events.txt", sep="\t", index_col=0, dtype=dtypes)

    # SNP scenario data
    d = {
        "chrom" : [],
        "motif_start" : [],
        "motif_end" : [],
        "strand" : [],
        "up_end" : [],
        "rtype" : [],
        "phase" : [],
        "splice_type" : [],
        "up_seq" : [],
        "ref_seq" : [],
        "alt_seq" : [],
        "down_seq" : [],
        "vnt_id" : [],
        "vnt_coord" : [],
        "vnt_pos" : [],
        "ref" : [],
        "alt" : [],
        "vnt_effect" : []
    }
    d.update({c : [] for c in maf_col[db]["dst"]})
    d.update({
        "vsite_ind" : [],
        "nsite_ind" : [],
        "ssite_ind" : [],
        "csite_inds" : [],
        "event_ind" : []
    })

    # for each SNP event:
    for e,event in events.iterrows():

        # cleavage sites
        cs = {int(p) : str(c) for p,c in zip(event["csite_pos"].split(","), event["csite_inds"].split(","))}

        # splice instances: (up_end, rtype, frame, up_seq)
        insts = {0 : [], 1 : []}
        insts.update({p : list(zip(csites.loc[c]["up_end"].split(","),
                                   csites.loc[c]["rtype"].split(","),
                                   csites.loc[c]["frame"].split(","),
                                   [s[-3:] for s in csites.loc[c]["scen_eflank_5ss"].split(",")]))
                      for p,c in cs.items()})

        # for all splice scenarios:
        for scen in set(insts[0] + insts[1]):
            
            # location information (first base after motif)
            [d[c].append(event[c]) for c in ["chrom", "motif_start", "motif_end", "strand"]]
            
            # scenario information
            d["up_end"].append(scen[0])
            d["rtype"].append(scen[1])
            d["phase"].append(scen[2])

            # alternatively-spliced
            if (scen in insts[0]) and (scen in insts[1]):
                d["splice_type"].append("AS")
                d["csite_inds"].append(f"{cs[0]},{cs[1]}")

            # proximally-spliced
            elif (scen in insts[0]):
                d["splice_type"].append("PS")
                d["csite_inds"].append(cs[0])
            
            # distally-spliced
            elif (scen in insts[1]):
                d["splice_type"].append("DS")
                d["csite_inds"].append(cs[1])

            # sequences
            d["up_seq"].append(scen[3])
            [d[c].append(event[c]) for c in ["ref_seq", "alt_seq", "down_seq"]]

            # SNP info: ID, coordinate, position, ref, alt, type
            [d[c].append(event[c]) for c in ["vnt_id", "vnt_coord", "vnt_pos", "ref", "alt", "vnt_effect"]]

            # minor allele frequencies
            [d[c].append(event[c]) for c in maf_col[db]["dst"]]

            # reference indices
            [d[c].append(event[c]) for c in ["vsite_ind", "nsite_ind", "ssite_ind"]]
            d["event_ind"].append(e)

    # DataFrame
    scens = pd.DataFrame(d)
    scens.to_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_scenarios.txt", sep="\t", index_label="index")

    print("\n===== SNP Scenarios =====")
    print(scens)

    # add scenario indices to SNPs
    snps["scen_inds"] = [",".join(map(str, scens[scens["vnt_id"] == i].index.values)) for i in snps.index]
    snps.to_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snps.txt", sep="\t")

    print("\n===== SNPs =====")
    print(snps)

    # add scenario indices to events
    events["scen_inds"] = [",".join(map(str, scens[scens["event_ind"] == i].index.values)) for i in events.index]
    events.to_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_events.txt", sep="\t")

    print("\n===== SNP Events =====")
    print(events)

''' Statistics '''

# Event + Scenario Frequencies
def count_events_scens():

    ''' Count NAGNAG-affecting SNP scenarios
    
    src:
    * events: /mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_events.txt
    * scenarios: /mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_scenarios.txt

    dst:
    * /mnt/data_3/hallep/nagnag/snps/statistics/nagnag_snp_events_scenarios_frequencies.txt
        * rows: (type, effect)
            * type: "events", "scenarios"
            * effect: "Create", "Alter", "Destroy", "All"
        * columns: {snp_dbs}
    * /mnt/data_3/hallep/nagnag/snps/statistics/nagnag_snp_scenarios_frequencies.txt
        * rows: (snp_db, effect)
            * snp_db: {snp_dbs}
            * effect: "Create", "Alter", "Destroy", "All"
        * columns: "Non-Coding", "Phase 0", "Phase 1", "Phase 2", "CDS", "All"
    
    Parameters
    ----------
    db : str {"dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice}
        SNP database
    '''

    # load SNP events + scenarios
    events = {db : pd.read_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_events.txt",
                               sep="\t", index_col=0, dtype=dtypes) for db in snp_dbs}
    scens = {db : pd.read_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_scenarios.txt",
                              sep="\t", index_col=0, dtype=dtypes) for db in snp_dbs}

    effects = ["CREATE", "ALTER", "DESTROY"]

    # events
    cE = np.array([[(df["vnt_effect"] == e).sum() for df in events.values()] for e in effects])
    cE = np.append(cE, np.sum(cE, axis=0).reshape(1, -1), axis=0)
    dfE = pd.DataFrame(cE, index=["Create", "Alter", "Destroy", "All"], columns=snp_dbs)

    # scenarios
    cS = np.array([[(df["vnt_effect"] == e).sum() for df in scens.values()] for e in effects])
    cS = np.append(cS, np.sum(cS, axis=0).reshape(1, -1), axis=0)
    dfS = pd.DataFrame(cS, index=["Create", "Alter", "Destroy", "All"], columns=snp_dbs)

    # combined
    df = pd.concat((dfE, dfS), axis=0, keys=["events", "scenarios"])
    df.to_csv("/mnt/data_3/hallep/nagnag/snps/statistics/nagnag_snp_events_scenarios_frequencies.txt",
              sep="\t", index_label=["type", "effect"])

    print("\n===== SNP Events + Scenarios =====")
    print(df)

    # by phase
    effects = ["CREATE", "ALTER", "DESTROY"]
    phases = [-1, 0, 1, 2]

    def by_phase(db:str) -> pd.DataFrame:

        df = scens[db]

        # count by phase + effect
        c = [[df[(df["phase"] == p) & (df["vnt_effect"] == e)].__len__() for p in phases] for e in effects]
        c = np.append(c, np.sum(c, axis=0).reshape(1, -1), axis=0)
        c = np.append(c, np.sum(c[:, 1:], axis=1).reshape(-1, 1), axis=1)
        c = np.append(c, np.sum(c[:, :4], axis=1).reshape(-1, 1), axis=1)

        # create DataFrame
        freq = pd.DataFrame(c, index=["Create", "Alter", "Destroy", "All"],
                            columns=["NC", "Phase 0", "Phase 1", "Phase 2", "CDS", "All"])

        return freq

    # by phase
    dfs = [by_phase(db) for db in snp_dbs]
    df = pd.concat(dfs, axis=0, keys=snp_dbs)
    df.to_csv("/mnt/data_3/hallep/nagnag/snps/statistics/nagnag_snp_scenarios_frequencies.txt",
              sep="\t", index_label=["snp_db", "effect"])

    print("\n===== SNP Scenarios =====")
    print(df)

# Affected Base Frequencies
def count_affected_bases():

    ''' Count NAGNAG-affecting SNP bases
    
    source: /mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_scenarios.txt

    destination: /mnt/data_3/hallep/nagnag/snps/statistics/affected_base.txt
        * rows: (effect, database)
            * effect: "CREATE", "ALTER", "DESTROY"
            * database: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice
        * columns: "N1", "A1", "G1", "N2", "A2", "G2"
    '''

    # load SNP events
    events = [pd.read_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_events.txt",
                          sep="\t", index_col=0, dtype=dtypes) for db in snp_dbs]

    # possible positions for each SNP effect type
    poss = {
        "CREATE" : [1, 2, 4, 5],
        "ALTER" : [0, 3],
        "DESTROY" : [1, 2, 4, 5]
    }

    # count frequencies
    counts = [[df[(df["vnt_effect"] == e) & (df["vnt_pos"] == p)].__len__()
               for df in events] for e,pos in poss.items() for p in pos]

    mpos = {0 : "N1", 1 : "A1", 2 : "G1", 3 : "N2", 4 : "A2", 5 : "G2"}
    [(e.title(), mpos[p]) for e,pos in poss.items() for p in pos]


    # DataFrame
    df = pd.DataFrame(counts, index=pd.MultiIndex.from_tuples([(e.title(), mpos[p])
                                                               for e,pos in poss.items() for p in pos],
                                                               names=["effect", "base"]),
                      columns=snp_dbs)

    df.to_csv("/mnt/data_3/hallep/nagnag/snps/statistics/affected_base_frequencies.txt", sep="\t")

    print("\n===== Affected Base =====")
    print(df)

''' ===== Proteome Effects ===== '''

# PROTEOME EFFECTS
def get_proteome_effects(db:str):

    ''' Get proteome effects for SNP scenarios

    scenarios: /mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_scenarios.txt
        * added columns: [metric]_[allele]
            * metric: "vc_ps", "vc_ds", "aa_ps", "aa_ds", "aattype",
                "ins_aa", "del_aa", "stop_lost", "stop_gained",
                "insert", "delete", "double", "hydrophobic_disrupt", "hydrophilic_disrupt",
            * allele: "ref", "alt" 

    Parameters
    ----------
    db : str {"dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice}
        SNP database
    '''

    # load scenarios
    scens = pd.read_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_scenarios.txt",
                        sep="\t", index_col=0, dtype=dtypes)

    # columns
    up = scens["up_seq"].values
    ref = scens["ref_seq"].values
    alt = scens["alt_seq"].values
    down = scens["down_seq"].values
    phase = scens["phase"].values

    # amino acid transitions
    tsn_ref = list(zip(*[get_transition(up=u, motif=m, down=d, phase=p, na0=False) for u,m,d,p in zip(up, ref, down, phase)]))
    tsn_alt = list(zip(*[get_transition(up=u, motif=m, down=d, phase=p, na0=False) for u,m,d,p in zip(up, alt, down, phase)]))

    for tr,ta,n in zip(tsn_ref, tsn_alt, TSNS):
        scens[f"{n}_ref"] = tr
        scens[f"{n}_alt"] = ta

    _, _, aa_ps_ref, aa_ds_ref, _ = [t for t in tsn_ref]
    _, _, aa_ps_alt, aa_ds_alt, _ = [t for t in tsn_alt]

    # binary amino acid outcomes
    binary_ref = zip(*[get_binary_outcomes(p, ps, ds) for p,ps,ds in zip(phase, aa_ps_ref, aa_ds_ref)])
    binary_alt = zip(*[get_binary_outcomes(p, ps, ds) for p,ps,ds in zip(phase, aa_ps_alt, aa_ds_alt)])

    for br,ba,n in zip(binary_ref, binary_alt, BIN):
        scens[f"{n}_ref"] = br
        scens[f"{n}_alt"] = ba

    # categorical amino acid outcomes
    cat_ref = zip(*[get_categorical_outcomes(p, ps, ds) for p,ps,ds in zip(phase, aa_ps_ref, aa_ds_ref)])
    cat_alt = zip(*[get_categorical_outcomes(p, ps, ds) for p,ps,ds in zip(phase, aa_ps_alt, aa_ds_alt)])
    
    for cr,ca,n in zip(cat_ref, cat_alt, CAT):
        scens[f"{n}_aa_ref"] = cr
        scens[f"{n}_aa_alt"] = ca

    # save DataFrame
    scens.to_csv(f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_scenarios.txt", sep="\t", index=True)

    print(scens)

''' Amino Acid Transitions '''

# VARIABLE CODON PROBABILITY
def vc_weight(phase:int, wsource:int|list[int]=0, stype:str="all"):

    ''' Get weighted variable codon probabilities

    Parameters
    ----------
    phase : int {0, 1, 2}
        NAGNAG phase
    wsource : list of ints (default = 0 --> all stochastic)
        weighting source
        * 0: stochastic (all 1/4)
        * 1: weighted
            * exonic basese from 1-NAGs
            * motif bases from NAGNAGs
    stype : str {"all", "AS", "PS", "DS"} (default = "all")
        splice type (for NAGNAGs) from which to get probabilities
    
        Returns
        -------
        _ : list[float]
            variable codon probabilities
    '''

    if isinstance(wsource, int):
        wsource = [wsource] if phase == 0 else [wsource] * 4
    
    # position-wise base probabilities
    p = {
        # Stochastic
        0 : pd.DataFrame(np.full((4, (5*4*6)), 1/4), index=N,
                         columns=pd.MultiIndex.from_product((["all", "nc", "p0", "p1", "p2"],
                                                             ["all", "AS", "PS", "DS"],
                                                             ["u-2", "u-1", "n1", "n2", "d1", "d2"]))),

        # weighted from 1-NAGs
        1 : pd.read_csv("/mnt/data_3/hallep/nagnag/splice_selection/1nag_poswise_prop.txt", sep="\t", index_col=0, header=[0,1,2]),

        # weighted from NAGNAGs        
        2 : pd.read_csv("/mnt/data_3/hallep/nagnag/splice_selection/nagnag_poswise_prop.txt", sep="\t", index_col=0, header=[0,1,2])
    }

    # phase
    phases = {
        0 : ["nc"],
        1 : ["p1", "nc", "p1", "p1"],
        2 : ["p2", "p2", "nc", "p2"]
    }

    # base positions
    pos = {
        0 : ["n2"],
        1 : ["u-1", "n2", "d1", "d2"],
        2 : ["u-2", "u-1", "n2", "d1"]
    }

    # weight sources
    wgts = {p : [0 if w == 0 else 2 if x == "n2" else 1 for w,x in zip(wsource, pos[p])] for p in range(3)}

    # splice types
    stypes = {p : [stype if pos == "n2" else "all" for pos in pos[p]] for p in range(3)}

    # variable bases
    variable = {
        0 : list(itertools.product(N)),
        1 : list(itertools.product(N, repeat=4)),
        2 : list(itertools.product(N, repeat=4))
    }

    return [np.prod([p[w][f][s][x][b] for w,f,s,x,b in zip(wgts[phase], phases[phase], stypes[phase], pos[phase], vb)]) for vb in variable[phase]]

# GET POSSIBLE VARIABLE CODONS
def get_possible_vc(phase_freq:tuple[int, int, int]=[1,1,1]) -> pd.DataFrame:

    ''' Get all possible NAGNAG variable codons

    src:
    * 1-NAGs: /mnt/data_3/hallep/nagnag/splice_selection/1nag_poswise_prop.txt
    * NAGNAGs: /mnt/data_3/hallep/nagnag/splice_selection/nagnag_poswise_prop.txt
    
    Parameters
    ----------
    phase_freq : tuple of 3 ints (default = [1, 1, 1])
        number of scenarios in each phase (for combining)
    
    Returns
    -------
    vc : pandas.DataFrame
        all possible variable codons in phases 0, 1, 2, and all
        * index: "phase", "vc_ps", "vc_ds"
        * columns: "aa_ps", "aa_ds", "aatype", "exp_stoch", "exp_splice", "exp_ps", "exp_ds", "exp_as", {BIN}, {CAT}
    '''

    # PHASE 0
    vc_ps, vc_ds, aa_ps, aa_ds, aatype = zip(*[get_transition(up="NNN", motif=f"NNN{n}AG", down="NNN", phase=0) for n in N])
    p0 = pd.DataFrame({
        "phase" : [0] * len(N),
        "vc_ps" : vc_ps,
        "vc_ds" : vc_ds,
        "aa_ps" : aa_ps,
        "aa_ds" : aa_ds,
        "aatype" : aatype,
        "exp_stoch" : vc_weight(phase=0, wsource=0, stype="all"),
        "exp_splice" : vc_weight(phase=0, wsource=1, stype="all"),
        "exp_ps" : vc_weight(phase=0, wsource=1, stype="PS"),
        "exp_ds" : vc_weight(phase=0, wsource=1, stype="DS"),
        "exp_as" : vc_weight(phase=0, wsource=1, stype="AS"),
    })

    # variable bases
    vb = list(itertools.product(N, repeat=4))

    # PHASE 1
    vc_ps, vc_ds, aa_ps, aa_ds, aatype = zip(*[get_transition(up=f"NN{v[0]}", motif=f"NNN{v[1]}AG", down=f"{v[2]}{v[3]}N", phase=1) for v in vb])
    p1 = pd.DataFrame({
        "phase" : [1] * len(vb),
        "vc_ps" : vc_ps,
        "vc_ds" : vc_ds,
        "aa_ps" : aa_ps,
        "aa_ds" : aa_ds,
        "aatype" : aatype,
        "exp_stoch" : vc_weight(phase=1, wsource=0, stype="all"),
        "exp_splice" : vc_weight(phase=1, wsource=[1, 1, 1, 1], stype="all"),
        "exp_ps" : vc_weight(phase=1, wsource=[1, 1, 1, 1], stype="PS"),
        "exp_ds" : vc_weight(phase=1, wsource=[1, 1, 1, 1], stype="DS"),
        "exp_as" : vc_weight(phase=1, wsource=[1, 1, 1, 1], stype="AS"),
    })

    # PHASE 2
    vc_ps, vc_ds, aa_ps, aa_ds, aatype = zip(*[get_transition(up=f"N{v[0]}{v[1]}", motif=f"NNN{v[2]}AG", down=f"{v[3]}NN", phase=2) for v in vb])
    p2 = pd.DataFrame({
        "phase" : [2] * len(vb),
        "vc_ps" : vc_ps,
        "vc_ds" : vc_ds,
        "aa_ps" : aa_ps,
        "aa_ds" : aa_ds,
        "aatype" : aatype,
        "exp_stoch" : vc_weight(phase=2, wsource=0, stype="all"),
        "exp_splice" : vc_weight(phase=2, wsource=[1, 1, 1, 1], stype="all"),
        "exp_ps" : vc_weight(phase=2, wsource=[1, 1, 1, 1], stype="PS"),
        "exp_ds" : vc_weight(phase=2, wsource=[1, 1, 1, 1], stype="DS"),
        "exp_as" : vc_weight(phase=2, wsource=[1, 1, 1, 1], stype="AS"),
    })

    # add amino acid outcomes
    for df in [p0, p1, p2]:

        # categorical outcomes
        categorical = zip(*[get_categorical_outcomes(p, ps, ds) for p,ps,ds in zip(df["phase"], df["aa_ps"], df["aa_ds"])])
        for n,c in zip(CAT, categorical):
            df[n] = c

        # binary outcomes
        binary = zip(*[get_binary_outcomes(p, ps, ds) for p,ps,ds in zip(df["phase"], df["aa_ps"], df["aa_ds"])])
        for n,b in zip(BIN, binary):
            df[n] = b

    # COMBINED
    p3 = pd.concat([p0, p1, p2])
    f0, f1, f2 = np.divide(phase_freq, sum(phase_freq))

    p3["phase"] = [3] * len(p3)
    for col in ["exp_stoch", "exp_splice", "exp_ps", "exp_ds", "exp_as"]:
        p3[col] = pd.concat([p0[col] * f0, p1[col] * f1, p2[col] * f2])

    # ALL
    vc = pd.concat([p0, p1, p2, p3]).set_index(keys=["phase", "vc_ps", "vc_ds"])

    return vc

''' Count '''

# TRANSITIONS
def count_transitions(db:str=None):

    ''' Find the frequencies of NAGNAG transitions
    
    src: 
    * canonical: /mnt/data_3/hallep/nagnag/nagnag_scenarios.txt
    * SNP: /mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_scenarios.txt

    EXP columns:
    * canonical: "exp_stoch", "exp_splice", "exp_as", "exp_ps", "exp_ds"
    * SNP: "exp_stoch", "exp_splice"

    FREQ columns:
    * canonical: "num", "num_as", "num_ps", "num_ds", "prop", "prop_as", "prop_ps", "prop_ds"
    * SNP: "num", "num_created", "num_altered_ref", "num_altered_alt", "num_destroyed",
           "prop", "prop_created", "prop_altered_ref", "prop_altered_alt", "prop_destroyed"

    variable codons
        * filepaths:
            * canonical: /mnt/data_3/hallep/nagnag/scenarios/nagnag_vc.txt
            * SNP: /mnt/data_3/hallep/nagnag/scenarios/{db}_vc.txt
        * index: "phase" (0-3), "vc_ps", "vc_ds"
        * columns: "aa_ps", "aa_ds", "aatype", {EXP}, {BIN}, {CAT}, {FREQ}

    amino acid transitions:
        * filepaths:
            * canonical: /mnt/data_3/hallep/nagnag/scenarios/nagnag_aat.txt
            * SNP: /mnt/data_3/hallep/nagnag/scenarios/{db}_aat.txt
        * index: "phase" (0-3), "aa_ps", "aa_ds"
        * columns: "num_vc", "prop_vc", "aatype", {EXP}, {BIN}, {CAT}, {FREQ}

    amino acid transition types:
        * filepaths:
            * canonical: /mnt/data_3/hallep/nagnag/scenarios/nagnag_aatt.txt
            * SNP: /mnt/data_3/hallep/nagnag/scenarios/{db}_aatt.txt
        * index: "phase" (0-2), "aatype"
        * columns: "num_vc", "prop_vc", "num_aat", "prop_aat", {EXP}, {BIN}, {CAT}, {FREQ}

    Parameters
    ----------
    db : str {"dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice} (default = None)
        SNP database \\
        if None, use canonical
    '''

    # define {src}, {vc_cols}, {dst_vc}, {dst_aat}, {dst_aatt}, {exp_col}, and {freq_col}
    if db:
        src = f"/mnt/data_3/hallep/nagnag/snps/{db}_nagnag_snp_scenarios.txt"
        vc_cols = ["vc_ps_ref", "vc_ps_alt", "vc_ds_ref", "vc_ds_alt"]
        
        dst_vc = f"/mnt/data_3/hallep/nagnag/scenarios/{db}_vc.txt"
        dst_aat = f"/mnt/data_3/hallep/nagnag/scenarios/{db}_aat.txt"
        dst_aatt = f"/mnt/data_3/hallep/nagnag/scenarios/{db}_aatt.txt"
        
        exp_col = ["exp_stoch", "exp_splice"]
        freq_col = ["num", "num_created", "num_altered_ref", "num_altered_alt", "num_destroyed",
                    "prop", "prop_created", "prop_altered_ref", "prop_altered_alt", "prop_destroyed"]
    else:
        src = "/mnt/data_3/hallep/nagnag/nagnag_scenarios.txt"
        vc_cols = ["vc_ps", "vc_ds"]
        
        dst_vc = "/mnt/data_3/hallep/nagnag/scenarios/nagnag_vc.txt"
        dst_aat = "/mnt/data_3/hallep/nagnag/scenarios/nagnag_aat.txt"
        dst_aatt = "/mnt/data_3/hallep/nagnag/scenarios/nagnag_aatt.txt"

        exp_col = ["exp_stoch", "exp_splice", "exp_as", "exp_ps", "exp_ds"]
        freq_col = ["num", "num_as", "num_ps", "num_ds", "prop", "prop_as", "prop_ps", "prop_ds"]

    # load scenarios
    scens = pd.read_csv(src, sep="\t", index_col=0, dtype=dtypes)
    scens = scens[scens["phase"] != -1]

    for c in vc_cols:
        scens[c] = scens[c].str.upper()

    scens_phase = [scens[scens["phase"] == p] for p in range(3)] + [scens]

    # get scenario variable codons
    if db:
        scen_ref = [list(zip(df["vnt_effect"], df["vc_ps_ref"], df["vc_ds_ref"])) for df in scens_phase]
        scen_alt = [list(zip(df["vnt_effect"], df["vc_ps_alt"], df["vc_ds_alt"])) for df in scens_phase]
    else:
        scen_vc = [list(zip(df["splice_type"], df["vc_ps"], df["vc_ds"])) for df in scens_phase]

    # VARIABLE CODONS
    vc = get_possible_vc(phase_freq=[len(scens[scens["phase"] == p]) for p in range(3)])
    codons = [vc.loc[p].index.values for p in range(4)]

    if db:
        vc.drop(columns=["exp_ps", "exp_ds", "exp_as"], inplace=True)

    # counts and proportions
    if db:
        # counts
        vc["num_created"] = [svcs.count(("CREATE",p,d)) for svcs,vcs in zip(scen_alt, codons) for p,d in vcs]
        vc["num_altered_ref"] = [svcs.count(("ALTER",p,d)) for svcs,vcs in zip(scen_ref, codons) for p,d in vcs]
        vc["num_altered_alt"] = [svcs.count(("ALTER",p,d)) for svcs,vcs in zip(scen_alt, codons) for p,d in vcs]
        vc["num_destroyed"] = [svcs.count(("DESTROY",p,d)) for svcs,vcs in zip(scen_ref, codons) for p,d in vcs]
        vc.insert(loc=len(vc.columns)-3, column="num", value=vc["num_created"] + vc["num_altered_ref"] + vc["num_altered_alt"] + vc["num_destroyed"])
        
        # proportions
        vc["prop"] = [x / vc.loc[p]["num"].sum() for p in range(4) for x in vc.loc[p]["num"]]
        vc["prop_created"] = [x / vc.loc[p]["num_created"].sum() for p in range(4) for x in vc.loc[p]["num_created"]]
        vc["prop_altered_ref"] = [x / vc.loc[p]["num_altered_ref"].sum() for p in range(4) for x in vc.loc[p]["num_altered_ref"]]
        vc["prop_altered_alt"] = [x / vc.loc[p]["num_altered_alt"].sum() for p in range(4) for x in vc.loc[p]["num_altered_alt"]]
        vc["prop_destroyed"] = [x / vc.loc[p]["num_destroyed"].sum() for p in range(4) for x in vc.loc[p]["num_destroyed"]]
    else:
        # counts
        vc["num_ps"] = [svcs.count(("PS",p,d)) for svcs,vcs in zip(scen_vc, codons) for p,d in vcs]
        vc["num_ds"] = [svcs.count(("DS",p,d)) for svcs,vcs in zip(scen_vc, codons) for p,d in vcs]
        vc["num_as"] = [svcs.count(("AS",p,d)) for svcs,vcs in zip(scen_vc, codons) for p,d in vcs]
        vc.insert(loc=len(vc.columns)-3, column="num", value=vc["num_ps"] + vc["num_ds"] + vc["num_as"])

        # proportions
        vc["prop"] = [x / vc.loc[p]["num"].sum() for p in range(4) for x in vc.loc[p]["num"]]
        vc["prop_as"] = [x / vc.loc[p]["num_as"].sum() for p in range(4) for x in vc.loc[p]["num_as"]]
        vc["prop_ps"] = [x / vc.loc[p]["num_ps"].sum() for p in range(4) for x in vc.loc[p]["num_ps"]]
        vc["prop_ds"] = [x / vc.loc[p]["num_ds"].sum() for p in range(4) for x in vc.loc[p]["num_ds"]]
    
    print("\n===== Variable Codons =====")
    print(vc)

    # save
    vc.to_csv(dst_vc, sep="\t")

    # AMINO ACID TRANSITIONS
    vc["phase"] = vc.index.get_level_values(0)
    vc.insert(loc=3, column="num_vc", value=1)
    vc.insert(loc=4, column="prop_vc", value=[1/len(vc.loc[p]) for p in range(4) for _ in range(len(vc.loc[p]))])
    vc_phase = [vc.loc[p] for p in range(4)]

    agg = {"aatype" : (lambda ser: ser.iloc[0])}
    agg.update({col : "sum" for col in ["num_vc", "prop_vc"] + exp_col})
    agg.update({col : (lambda ser: ser.iloc[0]) for col in BIN + CAT})
    agg.update({col : "sum" for col in freq_col})
    
    aat_phase = [v.groupby(["phase", "aa_ps", "aa_ds"]).aggregate(agg) for v in vc_phase]
    aat = pd.concat(aat_phase)

    print("\n===== Amino Acid Transitions =====")
    print(aat)

    # save
    aat.to_csv(dst_aat, sep="\t")

    # AMINO ACID TRANSITION TYPES
    aat["phase"] = aat.index.get_level_values(0)
    aat.insert(loc=3, column="num_aat", value=1)
    aat.insert(loc=4, column="prop_aat", value=[1/len(aat.loc[p]) for p in range(4) for _ in range(len(aat.loc[p]))])
    aat_phase = [aat.loc[p] for p in [0,1,2]]

    agg = {col : "sum" for col in ["num_vc", "prop_vc", "num_aat", "prop_aat"] + exp_col}
    agg.update({col : (lambda ser: ser.iloc[0]) for col in BIN + CAT})
    agg.update({col : "sum" for col in freq_col})
    
    aat0 = aat.loc[0].drop(columns=["aatype", "phase"])
    aat0.index = pd.MultiIndex.from_product([[0], aat0.index.get_level_values(0)], names=["phase", "aatype"])
    aatt_phase = [a.groupby(["phase", "aatype"]).aggregate(agg) for a in aat_phase]
    aatt = pd.concat([aat0] + aatt_phase).reindex([(0,i) for i in ["E", "Q", "K", "*"]] + list(itertools.product([1,2], ["DID", "NID", "CID", "IDR", "NC", "ET"])))

    print("\n===== Amino Acid Transition Types =====")
    print(aatt)

    # save
    aatt.to_csv(dst_aatt, sep="\t")

# OUTCOMES
def count_outcomes(db:str=None):

    ''' Find the frequencies of canonical amino acid outcomes
    
    src: 
    * canonical: /mnt/data_3/hallep/nagnag/scenarios/nagnag_vc.txt
    * SNP: /mnt/data_3/hallep/nagnag/scenarios/{db}_vc.txt

    dst:
        * filepaths:
            * canonical: /mnt/data_3/hallep/nagnag/scenarios/nagnag_outcomes.txt
            * SNP: /mnt/data_3/hallep/nagnag/scenarios/{db}_outcomes.txt
        * index: phase, cat, val
            * phase: {0, 1, 2, 3}
            * cat: {"bin", "ins", "del"}
            * val: {BIN}, {AA}, {AA}
        * columns:
            * canonical: "exp_stoch", "exp_splice", "exp_as", "exp_ps", "exp_ds",
                    "num", "num_as", "num_ps", "num_ds", "prop", "prop_as", "prop_ps", "prop_ds"
            * SNP: "exp_stoch", "exp_splice", 
                    "num", "num_created", "num_altered_ref", "num_altered_alt", "num_destroyed",
                    "prop", "prop_created", "prop_altered_ref", "prop_altered_alt", "prop_destroyed"

    Parameters
    ----------
    db : str {"dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice} (default = None)
        SNP database \\
        if None, use canonical
    '''

    # define {src}, {dst}, and {columns}
    if db:
        src = f"/mnt/data_3/hallep/nagnag/scenarios/{db}_vc.txt"
        dst = f"/mnt/data_3/hallep/nagnag/scenarios/{db}_outcomes.txt"

        columns = ["exp_stoch", "exp_splice", "num", "num_created", "num_altered_ref", "num_altered_alt", "num_destroyed",
                   "prop", "prop_created", "prop_altered_ref", "prop_altered_alt", "prop_destroyed"]
    else:
        src = "/mnt/data_3/hallep/nagnag/scenarios/nagnag_vc.txt"
        dst = "/mnt/data_3/hallep/nagnag/scenarios/nagnag_outcomes.txt"

        columns = ["exp_stoch", "exp_splice", "exp_as", "exp_ps", "exp_ds",
                   "num", "num_as", "num_ps", "num_ds", "prop", "prop_as", "prop_ps", "prop_ds"]

    # load variable codons
    vc = pd.read_csv(src, sep="\t")
    ph = vc["phase"]

    # inserted/deleted amino acids
    cats = list(itertools.product(CAT, AA))

    # count outcomes
    cBin = {p : {col : [vc[(ph == p) & (vc[b] == 1)][col].sum() for b in BIN] for col in columns} for p in range(4)}
    cCat = {p : {col : [vc[(ph == p) & (vc[c].str.contains(a, regex=False))][col].sum()
                        for c,a in cats] for col in columns} for p in range(4)}

    # DataFrame
    index = [(a,b,c) for a,(b,c) in itertools.product(range(4), list(itertools.product(["bin"], BIN)) + cats)]
    df = pd.DataFrame({col : [x for p in range(4) for o in [cBin[p][col], cCat[p][col]] for x in o] for col in columns},
                      index=pd.MultiIndex.from_tuples(index, names=["phase", "cat", "val"]))

    print("\n===== Outcomes =====")
    print(df)

    df.to_csv(dst, sep="\t")

''' Aggregate '''

# VARIABLE CODONS
def aggregate_codons():

    ''' Combine all canonical and SNP variable codons
    
    src: 
    * canonical: /mnt/data_3/hallep/nagnag/scenarios/nagnag_vc.txt
    * SNP: /mnt/data_3/hallep/nagnag/scenarios/{db}_vc.txt

    proportions: /mnt/data_3/hallep/nagnag/scenarios/variable_codons.txt
        * index: (phase, vc_ps, vc_ds)
            * phase: {0, 1, 2}
        * columns: 
            * ("canon", freq)
                * freq: "exp_stoch", "exp_splice", "obs"
            * (splice type, freq)
                * splice type: "ps", "ds", "as"
                * freq: "exp", "obs"
                * e: "exp_stoch", "exp_splice", "exp_as", "exp_ps", "exp_ds"
            * (effect, db)
                * effect: "create", "alter_ref", "alter_alt", "destroy"
                * db: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"
        * sorted by phase (0, 1, 2) and canonical frequency (descending)
    
    counts: /mnt/data_3/hallep/nagnag/scenarios/variable_codons_counts.txt
        * index: (phase, aattype)
            * phase: {0, 1, 2}
        * columns: 
            * ("canon", splice type)
                * splice type: "all", "ps", "ds", "as"
            * (effect, db)
                * effect: "create", "alter_ref", "alter_alt", "destroy"
                * db: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"
        * sorted by phase (0, 1, 2) and canonical frequency (descending)
    '''

    # load transition types
    canon = pd.read_csv("/mnt/data_3/hallep/nagnag/scenarios/nagnag_vc.txt", sep="\t", index_col=[0,1,2]).loc[[0,1,2]]    
    dbs = ["dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"]
    snps = ({db : pd.read_csv(f"/mnt/data_3/hallep/nagnag/scenarios/{db}_vc.txt", sep="\t", index_col=[0,1,2]).loc[[0,1,2]] for db in dbs})

    # column name dictionary
    cd = {
        "exp_stoch" : "exp_stoch",
        "exp_splice" : "exp_splice",
        "exp" : "exp",
        "obs" : "prop",
        "all" : "",
        "ps" : "_ps",
        "ds" : "_ds",
        "as" : "_as",
        "create" : "created",
        "alter_ref" : "altered_ref",
        "alter_alt" : "altered_alt",
        "destroy" : "destroyed"
    }

    # Proportions
    p = {("canon", f) : canon[cd[f]].values for f in ["exp_stoch", "exp_splice", "obs"]}
    p.update({(st, f) : canon[f"{cd[f]}_{st}"].values for st,f in itertools.product(["ps", "ds", "as"], ["exp", "obs"])})
    p.update({(e, db) : df[f"prop_{cd[e]}"].values for e,(db,df) in itertools.product(nn_effects, snps.items())})

    # create DataFrame
    dfP = pd.DataFrame(p, index=canon.index)
    dfP[("i","j")] = dfP.index.get_level_values(0).values
    dfP.sort_values(by=[("i", "j"), ("canon", "obs")], ascending=[True, False], inplace=True)
    dfP.drop(columns=("i","j"), inplace=True)
    dfP.to_csv("/mnt/data_3/hallep/nagnag/scenarios/variable_codons.txt", sep="\t", index=True)

    print("\n===== Amino Acid Transitions =====")
    print(dfP)

    # Counts
    c = {("canon", s) : canon[f"num{cd[s]}"] for s in ["all", "ps", "ds", "as"]}
    c.update({(e, db) : df[f"num_{cd[e]}"].values for e,(db,df) in itertools.product(nn_effects, snps.items())})

    dfC = pd.DataFrame(c, index=canon.index)
    dfC[("i","j")] = dfC.index.get_level_values(0).values
    dfC.sort_values(by=[("i", "j"), ("canon", "all")], ascending=[True, False], inplace=True)
    dfC.drop(columns=("i","j"), inplace=True)
    dfC.to_csv("/mnt/data_3/hallep/nagnag/scenarios/variable_codons_counts.txt", sep="\t")

    print("\n===== Amino Acid Transitions (Raw Counts) =====")
    print(dfC)

# TRANSITIONS
def aggregate_transitions():

    ''' Combine all canonical and SNP amino acid transitions
    
    src: 
    * canonical: /mnt/data_3/hallep/nagnag/scenarios/nagnag_aat.txt
    * SNP: /mnt/data_3/hallep/nagnag/scenarios/{db}_aat.txt

    proportions: /mnt/data_3/hallep/nagnag/scenarios/aa_transitions.txt
        * index: (phase, aa_ps, aa_ds)
            * phase: {0, 1, 2}
        * columns: 
            * ("canon", freq)
                * freq: "exp_stoch", "exp_splice", "obs"
            * (splice type, freq)
                * splice type: "ps", "ds", "as"
                * freq: "exp", "obs"
                * e: "exp_stoch", "exp_splice", "exp_as", "exp_ps", "exp_ds"
            * (effect, db)
                * effect: "create", "alter_ref", "alter_alt", "destroy"
                * db: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"
        * sorted by phase (0, 1, 2) and canonical frequency (descending)
    
    counts: /mnt/data_3/hallep/nagnag/scenarios/aa_transitions_counts.txt
        * index: (phase, aattype)
            * phase: {0, 1, 2}
        * columns: 
            * ("canon", splice type)
                * splice type: "all", "ps", "ds", "as"
            * (effect, db)
                * effect: "create", "alter_ref", "alter_alt", "destroy"
                * db: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"
        * sorted by phase (0, 1, 2) and canonical frequency (descending)
    '''

    # load transition types
    canon = pd.read_csv("/mnt/data_3/hallep/nagnag/scenarios/nagnag_aat.txt", sep="\t", index_col=[0,1,2]).loc[[0,1,2]]    
    dbs = ["dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"]
    snps = ({db : pd.read_csv(f"/mnt/data_3/hallep/nagnag/scenarios/{db}_aat.txt", sep="\t", index_col=[0,1,2]).loc[[0,1,2]] for db in dbs})

    # column name dictionary
    cd = {
        "exp_stoch" : "exp_stoch",
        "exp_splice" : "exp_splice",
        "exp" : "exp",
        "obs" : "prop",
        "all" : "",
        "ps" : "_ps",
        "ds" : "_ds",
        "as" : "_as",
        "create" : "created",
        "alter_ref" : "altered_ref",
        "alter_alt" : "altered_alt",
        "destroy" : "destroyed"
    }

    # Proportions
    p = {("canon", f) : canon[cd[f]].values for f in ["exp_stoch", "exp_splice", "obs"]}
    p.update({(st, f) : canon[f"{cd[f]}_{st}"].values for st,f in itertools.product(["ps", "ds", "as"], ["exp", "obs"])})
    p.update({(e, db) : df[f"prop_{cd[e]}"].values for e,(db,df) in itertools.product(nn_effects, snps.items())})

    # create DataFrame
    dfP = pd.DataFrame(p, index=canon.index)
    dfP[("i","j")] = dfP.index.get_level_values(0).values
    dfP.sort_values(by=[("i", "j"), ("canon", "obs")], ascending=[True, False], inplace=True)
    dfP.drop(columns=("i","j"), inplace=True)
    dfP.to_csv("/mnt/data_3/hallep/nagnag/scenarios/aa_transitions.txt", sep="\t", index=True)

    print("\n===== Amino Acid Transitions =====")
    print(dfP)

    # Counts
    c = {("canon", s) : canon[f"num{cd[s]}"] for s in ["all", "ps", "ds", "as"]}
    c.update({(e, db) : df[f"num_{cd[e]}"].values for e,(db,df) in itertools.product(nn_effects, snps.items())})

    dfC = pd.DataFrame(c, index=canon.index)
    dfC[("i","j")] = dfC.index.get_level_values(0).values
    dfC.sort_values(by=[("i", "j"), ("canon", "all")], ascending=[True, False], inplace=True)
    dfC.drop(columns=("i","j"), inplace=True)
    dfC.to_csv("/mnt/data_3/hallep/nagnag/scenarios/aa_transitions_counts.txt", sep="\t")

    print("\n===== Amino Acid Transitions (Raw Counts) =====")
    print(dfC)

# TRANSITION TYPES
def aggregate_transition_types():

    ''' Combine all canonical and SNP amino acid transition types
    
    src: 
    * canonical: /mnt/data_3/hallep/nagnag/scenarios/nagnag_aatt.txt
    * SNP: /mnt/data_3/hallep/nagnag/scenarios/{db}_aatt.txt

    proportions: /mnt/data_3/hallep/nagnag/scenarios/aa_transition_types.txt
        * index: (phase, aattype)
            * phase: {0, 1, 2}
        * columns: 
            * ("canon", freq)
                * freq: "exp_stoch", "exp_splice", "obs"
            * (splice type, freq)
                * splice type: "ps", "ds", "as"
                * freq: "exp", "obs"
                * e: "exp_stoch", "exp_splice", "exp_as", "exp_ps", "exp_ds"
            * (effect, db)
                * effect: "create", "alter_ref", "alter_alt", "destroy"
                * db: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"
    
    counts: /mnt/data_3/hallep/nagnag/scenarios/aa_transition_types_counts.txt
        * index: (phase, aattype)
            * phase: {0, 1, 2}
        * columns: 
            * ("canon", splice type)
                * splice type: "all", "ps", "ds", "as"
            * (effect, db)
                * effect: "create", "alter_ref", "alter_alt", "destroy"
                * db: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"
    '''

    # load transition types
    canon = pd.read_csv("/mnt/data_3/hallep/nagnag/scenarios/nagnag_aatt.txt", sep="\t", index_col=[0,1])
    dbs = ["dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"]
    snps = ({db : pd.read_csv(f"/mnt/data_3/hallep/nagnag/scenarios/{db}_aatt.txt", sep="\t", index_col=[0,1]) for db in dbs})

    # column name dictionary
    cd = {
        "exp_stoch" : "exp_stoch",
        "exp_splice" : "exp_splice",
        "exp" : "exp",
        "obs" : "prop",
        "all" : "",
        "ps" : "_ps",
        "ds" : "_ds",
        "as" : "_as",
        "create" : "created",
        "alter_ref" : "altered_ref",
        "alter_alt" : "altered_alt",
        "destroy" : "destroyed"
    }

    # Proportions
    p = {("canon", f) : canon[cd[f]].values for f in ["exp_stoch", "exp_splice", "obs"]}
    p.update({(st, f) : canon[f"{cd[f]}_{st}"].values for st,f in itertools.product(["ps", "ds", "as"], ["exp", "obs"])})
    p.update({(e, db) : df[f"prop_{cd[e]}"].values for e,(db,df) in itertools.product(nn_effects, snps.items())})

    # create DataFrame
    dfP = pd.DataFrame(p, index=canon.index)
    dfP.to_csv("/mnt/data_3/hallep/nagnag/scenarios/aa_transition_types.txt", sep="\t", index=True)

    print("\n===== Amino Acid Transitions =====")
    print(dfP)

    # Counts
    c = {("canon", s) : canon[f"num{cd[s]}"] for s in ["all", "ps", "ds", "as"]}
    c.update({(e, db) : df[f"num_{cd[e]}"].values for e,(db,df) in itertools.product(nn_effects, snps.items())})

    dfC = pd.DataFrame(c, index=canon.index)
    dfC.to_csv("/mnt/data_3/hallep/nagnag/scenarios/aa_transition_types_counts.txt", sep="\t")

    print("\n===== Amino Acid Transitions (Raw Counts) =====")
    print(dfC)

# OUTCOMES
def aggregate_outcomes():

    ''' Combine all canonical and SNP amino acid outcomes
    
    src:
    * canonical: /mnt/data_3/hallep/nagnag/scenarios/nagnag_outcomes.txt
    * SNP: /mnt/data_3/hallep/nagnag/scenarios/{db}_outcomes.txt

    proportions: /mnt/data_3/hallep/nagnag/scenarios/aa_outcomes{_by_phase}.txt; 
        * index: (cat, val)
            * cat: "bin", "ins", "del"
            * val: {BIN}, {AA}, {AA}
        * columns:
            * ("canon", freq)
                * freq: "exp_stoch", "exp_splice", "obs"
            * (splice type, freq)
                * splice type: "ps", "ds", "as"
                * freq: "exp", "obs"
                * e: "exp_stoch", "exp_splice", "exp_as", "exp_ps", "exp_ds"
            * (effect, db)
                * effect: "create", "alter_ref", "alter_alt", "destroy"
                * db: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"
    
    counts: /mnt/data_3/hallep/nagnag/scenarios/aa_outcomes_counts{_by_phase}.txt
        * index: (cat, val)
            * cat: "bin", "ins", "del"
            * val: {BIN}, {AA}, {AA}
        * columns: 
            * ("canon", splice type)
                * splice type: "all", "ps", "ds", "as"
            * (effect, db)
                * effect: "create", "alter_ref", "alter_alt", "destroy"
                * db: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"
    '''

    # load outcomes
    canon = pd.read_csv("/mnt/data_3/hallep/nagnag/scenarios/nagnag_outcomes.txt", sep="\t", index_col=[0,1,2])
    snps = {db : pd.read_csv(f"/mnt/data_3/hallep/nagnag/scenarios/{db}_outcomes.txt",
                             sep="\t", index_col=[0,1,2]) for db in snp_dbs}

    # column name dictionary
    cd = {
        "exp_stoch" : "exp_stoch",
        "exp_splice" : "exp_splice",
        "exp" : "exp",
        "obs" : "prop",
        "all" : "",
        "ps" : "_ps",
        "ds" : "_ds",
        "as" : "_as",
        "create" : "created",
        "alter_ref" : "altered_ref",
        "alter_alt" : "altered_alt",
        "destroy" : "destroyed"
    }

    # Proportions
    p = {("canon", f) : canon[cd[f]].values for f in ["exp_stoch", "exp_splice", "obs"]}
    p.update({(st, f) : canon[f"{cd[f]}_{st}"].values for st,f in itertools.product(["ps", "ds", "as"], ["exp", "obs"])})
    p.update({(e, db) : df[f"prop_{cd[e]}"].values for e,(db,df) in itertools.product(nn_effects, snps.items())})

    # DataFrame
    dfP = pd.DataFrame(p, index=canon.index)
    dfP.to_csv("/mnt/data_3/hallep/nagnag/scenarios/aa_outcomes_by_phase.txt", sep="\t")
    dfP.loc[3].to_csv("/mnt/data_3/hallep/nagnag/scenarios/aa_outcomes.txt", sep="\t")

    print("\n===== Amino Acid Outcomes =====")
    print(dfP)

    # Counts
    c = {("canon", s) : canon[f"num{cd[s]}"] for s in ["all", "ps", "ds", "as"]}
    c.update({(e, db) : df[f"num_{cd[e]}"].values for e,(db,df) in itertools.product(nn_effects, snps.items())})

    dfC = pd.DataFrame(c, index=canon.index)
    dfC.to_csv("/mnt/data_3/hallep/nagnag/scenarios/aa_outcomes_counts_by_phase.txt", sep="\t")
    dfC.loc[3].to_csv("/mnt/data_3/hallep/nagnag/scenarios/aa_outcomes_counts.txt", sep="\t")

    print("\n===== Amino Acid Outcomes (Raw Counts) =====")
    print(dfC)

''' ===== RUN ===== '''

def analyze_db(db:str):
    get_snp_events(db)
    get_snp_scenarios(db)

''' Observed '''
[get_snp_events(db) for db in snp_dbs]
[get_snp_scenarios(db) for db in snp_dbs]

''' Statistics '''
count_events_scens()
count_affected_bases()

''' Proteome '''
[get_proteome_effects(db) for db in snp_dbs]
[count_outcomes(db) for db in [None] + snp_dbs]
[count_transitions(db) for db in [None] + snp_dbs]

''' Aggregate '''
aggregate_codons()
aggregate_transitions()
aggregate_transition_types()
aggregate_outcomes()
