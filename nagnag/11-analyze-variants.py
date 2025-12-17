''' Analyze the effects of variants on the proteome '''

from lib import argparse, itertools, pd, np
from utils import ROOT, log_script, log_fn
from utils.vnt import dtypes, all_vnt_dbs
from utils.seq import N, TSNS, OUT, AA, get_tsn, get_categorical_outcomes

# Variant Databases
parser = argparse.ArgumentParser()
parser.add_argument("-H", "--ignore-HGMD", action="store_true", help="Do not process/analyze HGMD Splice data")
args = parser.parse_args()

vnt_dbs = all_vnt_dbs
if args.ignore_HGMD:
    vnt_dbs.remove("HGMD_splice")

# --- Proteome Effects --- #

# Get proteome effects
def get_proteome_effects(db:str):

    ''' Get proteome effects for variant scenarios

    `variants/found/{db}_nagnag_vnt_scenarios.txt`
    ----------------------------------------------
    added columns: "[metric]_[allele]"
    * metric: "vc_ps", "vc_ds", "aa_ps", "aa_ds", "aattype", "ins_aa", "del_aa"
    * allele: "ref", "alt" 
    '''

    # load scenarios
    scens = pd.read_csv(f"{ROOT}/variants/found/{db}_nagnag_vnt_scenarios.txt",
                        sep="\t", index_col=0, dtype=dtypes)

    # columns
    up = scens["up_seq"].values
    ref = scens["ref_seq"].values
    alt = scens["alt_seq"].values
    down = scens["down_seq"].values
    phase = scens["phase"].values

    # amino acid transitions
    tsn_ref = list(zip(*[get_tsn(up=u, motif=m, down=d, phase=p) for u,m,d,p in zip(up, ref, down, phase)]))
    tsn_alt = list(zip(*[get_tsn(up=u, motif=m, down=d, phase=p) for u,m,d,p in zip(up, alt, down, phase)]))

    for tr,ta,n in zip(tsn_ref, tsn_alt, TSNS):
        scens[f"{n}_ref"] = tr
        scens[f"{n}_alt"] = ta

    _, _, aa_ps_ref, aa_ds_ref, _ = [t for t in tsn_ref]
    _, _, aa_ps_alt, aa_ds_alt, _ = [t for t in tsn_alt]

    # amino acid outcomes
    cat_ref = zip(*[get_categorical_outcomes(p, ps, ds) for p,ps,ds in zip(phase, aa_ps_ref, aa_ds_ref)])
    cat_alt = zip(*[get_categorical_outcomes(p, ps, ds) for p,ps,ds in zip(phase, aa_ps_alt, aa_ds_alt)])
    
    for cr,ca,n in zip(cat_ref, cat_alt, OUT):
        scens[f"{n}_aa_ref"] = cr
        scens[f"{n}_aa_alt"] = ca

    # save DataFrame
    scens.to_csv(f"{ROOT}/variants/found/{db}_nagnag_vnt_scenarios.txt", sep="\t", index=True)

# --- Proteome Effect Frequencies --- #

# Variable codon expected frequency
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
        1 : pd.read_csv(f"{ROOT}/proteome/1nag_poswise_freq.txt", sep="\t", index_col=0, header=[0,1,2]),

        # weighted from NAGNAGs        
        2 : pd.read_csv(f"{ROOT}/proteome/nagnag_poswise_freq.txt", sep="\t", index_col=0, header=[0,1,2])
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

# Possible variable codons
def get_possible_vc(phase_freq:dict[int,dict[str,int]]) -> pd.DataFrame:

    ''' Get all possible NAGNAG variable codons

    src:
    * 1-NAGs: `proteome/1nag_poswise_freq.txt`
    * NAGNAGs: `proteome/nagnag_poswise_freq.txt`
    
    Parameters
    ----------
    phase_freq : dictionary of dictionaries of integers
        number of scenarios in each phase (outer dict.) and splice type (inner dict.)
    
    Returns
    -------
    vc : pandas.DataFrame
        all possible variable codons in phases 0, 1, 2, and all
        * index: "phase", "vc_ps", "vc_ds"
        * columns: "aa_ps", "aa_ds", "aattype", "exp_stoch", "exp_splice", "exp_ps", "exp_ds", "exp_as", {BIN}, {CAT}
    '''

    # PHASE 0
    vc_ps, vc_ds, aa_ps, aa_ds, aattype = zip(*[get_tsn(up="NNN", motif=f"NNN{n}AG", down="NNN", phase=0) for n in N])
    p0 = pd.DataFrame({
        "phase" : [0] * len(N),
        "vc_ps" : vc_ps,
        "vc_ds" : vc_ds,
        "aa_ps" : aa_ps,
        "aa_ds" : aa_ds,
        "aattype" : aattype,
        "exp_stoch" : vc_weight(phase=0, wsource=0, stype="all"),
        "exp_splice" : vc_weight(phase=0, wsource=1, stype="all"),
        "exp_ps" : vc_weight(phase=0, wsource=1, stype="PS"),
        "exp_ds" : vc_weight(phase=0, wsource=1, stype="DS"),
        "exp_as" : vc_weight(phase=0, wsource=1, stype="AS"),
    })

    # variable bases
    vb = list(itertools.product(N, repeat=4))

    # PHASE 1
    vc_ps, vc_ds, aa_ps, aa_ds, aattype = zip(*[get_tsn(up=f"NN{v[0]}", motif=f"NNN{v[1]}AG", down=f"{v[2]}{v[3]}N", phase=1) for v in vb])
    p1 = pd.DataFrame({
        "phase" : [1] * len(vb),
        "vc_ps" : vc_ps,
        "vc_ds" : vc_ds,
        "aa_ps" : aa_ps,
        "aa_ds" : aa_ds,
        "aattype" : aattype,
        "exp_stoch" : vc_weight(phase=1, wsource=0, stype="all"),
        "exp_splice" : vc_weight(phase=1, wsource=[1, 1, 1, 1], stype="all"),
        "exp_ps" : vc_weight(phase=1, wsource=[1, 1, 1, 1], stype="PS"),
        "exp_ds" : vc_weight(phase=1, wsource=[1, 1, 1, 1], stype="DS"),
        "exp_as" : vc_weight(phase=1, wsource=[1, 1, 1, 1], stype="AS"),
    })

    # PHASE 2
    vc_ps, vc_ds, aa_ps, aa_ds, aattype = zip(*[get_tsn(up=f"N{v[0]}{v[1]}", motif=f"NNN{v[2]}AG", down=f"{v[3]}NN", phase=2) for v in vb])
    p2 = pd.DataFrame({
        "phase" : [2] * len(vb),
        "vc_ps" : vc_ps,
        "vc_ds" : vc_ds,
        "aa_ps" : aa_ps,
        "aa_ds" : aa_ds,
        "aattype" : aattype,
        "exp_stoch" : vc_weight(phase=2, wsource=0, stype="all"),
        "exp_splice" : vc_weight(phase=2, wsource=[1, 1, 1, 1], stype="all"),
        "exp_ps" : vc_weight(phase=2, wsource=[1, 1, 1, 1], stype="PS"),
        "exp_ds" : vc_weight(phase=2, wsource=[1, 1, 1, 1], stype="DS"),
        "exp_as" : vc_weight(phase=2, wsource=[1, 1, 1, 1], stype="AS"),
    })

    # add amino acid outcomes
    for df in [p0, p1, p2]:
        outcomes = zip(*[get_categorical_outcomes(p, ps, ds) for p,ps,ds in zip(df["phase"], df["aa_ps"], df["aa_ds"])])
        for i,(n,c) in enumerate(zip(OUT, outcomes)):
            df.insert(loc=i+6, column=n, value=c)

    # COMBINED
    p3 = pd.concat([p0, p1, p2])

    f = np.array([list(ph.values()) for ph in phase_freq.values()])
    f0, f1, f2 = f.sum(axis=1) / f.sum()
    
    p3["phase"] = [3] * len(p3)
    for col in ["exp_stoch", "exp_splice"]:
        p3[col] = pd.concat([p0[col] * f0, p1[col] * f1, p2[col] * f2])

    a = f / f.sum(axis=1)
    for r,s in enumerate(["ps", "ds", "as"]):
        p3[f"exp_{s}"] = pd.concat([p0[col] * a[0][r], p1[col] * a[1][r], p2[col] * a[2][r]])

    # ALL
    vc = pd.concat([p0, p1, p2, p3]).set_index(keys=["phase", "vc_ps", "vc_ds"])

    return vc

# Count transitions
def count_transitions(db:str=None):

    ''' Find the frequencies of NAGNAG transitions
    
    **Source:** 
    * reference: `sites/nagnag_scens.txt`
    * variant: `variants/found/{db}_nagnag_vnt_scenarios.txt`

    EXP columns:
    * reference: "exp_stoch", "exp_splice", "exp_as", "exp_ps", "exp_ds"
    * variant: "exp_stoch", "exp_splice"

    FREQ columns:
    * reference: "num", "num_as", "num_ps", "num_ds", "prop", "prop_as", "prop_ps", "prop_ds"
    * variant: "num", "num_create", "num_alter_ref", "num_alter_alt", "num_destroy",
           "prop", "prop_create", "prop_alter_ref", "prop_alter_alt", "prop_destroy"

    `variants/proteome/{nagnag/db}_vc.txt`
    --------------------------------------
    variable codons
    * index: "phase" (0-3), "vc_ps", "vc_ds"
    * columns: "aa_ps", "aa_ds", "aattype", "ins_aa", "del_aa", EXP, FREQ

    `variants/proteome/{nagnag/db}_aat.txt`
    ---------------------------------------
    amino acid transitions
    * index: "phase" (0-3), "aa_ps", "aa_ds"
    * columns: "num_vc", "prop_vc", "aattype", "ins_aa", "del_aa", EXP, FREQ

    `variants/proteome/{nagnag/db}_aatt.txt`
    ----------------------------------------
    amino acid transition types
    * index: "phase" (0-2), "aattype"
    * columns: "num_vc", "prop_vc", "num_aat", "prop_aat", "ins_aa", "del_aa", EXP, FREQ
    '''

    # define {src}, {vc_cols}, {dst_vc}, {dst_aat}, {dst_aatt}, {exp_col}, and {freq_col}
    if db:
        src = f"{ROOT}/variants/found/{db}_nagnag_vnt_scenarios.txt"
        vc_cols = ["vc_ps_ref", "vc_ps_alt", "vc_ds_ref", "vc_ds_alt"]
        
        dst_vc = f"{ROOT}/variants/proteome/{db}_vc.txt"
        dst_aat = f"{ROOT}/variants/proteome/{db}_aat.txt"
        dst_aatt = f"{ROOT}/variants/proteome/{db}_aatt.txt"
        
        exp_col = ["exp_stoch", "exp_splice"]
        freq_col = ["num", "num_create", "num_alter_ref", "num_alter_alt", "num_destroy",
                    "prop", "prop_create", "prop_alter_ref", "prop_alter_alt", "prop_destroy"]
    else:
        src = f"{ROOT}/sites/nagnag_scens.txt"
        vc_cols = ["vc_ps", "vc_ds"]
        
        dst_vc = f"{ROOT}/variants/proteome/nagnag_vc.txt"
        dst_aat = f"{ROOT}/variants/proteome/nagnag_aat.txt"
        dst_aatt = f"{ROOT}/variants/proteome/nagnag_aatt.txt"

        exp_col = ["exp_stoch", "exp_splice", "exp_as", "exp_ps", "exp_ds"]
        freq_col = ["num", "num_as", "num_ps", "num_ds", "prop", "prop_as", "prop_ps", "prop_ds"]

    # load scenarios
    scens = pd.read_csv(src, sep="\t", index_col=0, dtype=dtypes)
    scens = scens[scens["phase"] != -1]

    for c in vc_cols:
        scens[c] = scens[c].str.upper()

    scens_phase = [scens[scens["phase"] == p] for p in range(3)] + [scens]

    # get observed scenario variable codons by phase
    if db:
        scen_ref = [list(zip(df["vnt_effect"], df["vc_ps_ref"], df["vc_ds_ref"])) for df in scens_phase]
        scen_alt = [list(zip(df["vnt_effect"], df["vc_ps_alt"], df["vc_ds_alt"])) for df in scens_phase]
    else:
        scen_vc = [list(zip(df["splice_type"], df["vc_ps"], df["vc_ds"])) for df in scens_phase]

    # Variable Codons
    vc = get_possible_vc({p : {s : len(scens[(scens["phase"] == p) & (scens["splice_type"] == s)]) 
                               for s in ["PS", "DS", "AS"]} for p in range(3)})
    
    # possible variable codons
    codons = [vc.loc[p].index.values for p in range(4)]

    if db:
        vc.drop(columns=["exp_ps", "exp_ds", "exp_as"], inplace=True)

    # counts and proportions
    if db:
        # counts
        vc["num_create"] = [svcs.count(("CREATE",p,d)) for svcs,vcs in zip(scen_alt, codons) for p,d in vcs]
        vc["num_alter_ref"] = [svcs.count(("ALTER",p,d)) for svcs,vcs in zip(scen_ref, codons) for p,d in vcs]
        vc["num_alter_alt"] = [svcs.count(("ALTER",p,d)) for svcs,vcs in zip(scen_alt, codons) for p,d in vcs]
        vc["num_destroy"] = [svcs.count(("DESTROY",p,d)) for svcs,vcs in zip(scen_ref, codons) for p,d in vcs]
        vc.insert(loc=len(vc.columns)-3, column="num", value=vc["num_create"] + vc["num_alter_ref"] + vc["num_alter_alt"] + vc["num_destroy"])
        
        # proportions
        vc["prop"] = [x / vc.loc[p]["num"].sum() for p in range(4) for x in vc.loc[p]["num"]]
        vc["prop_create"] = [x / vc.loc[p]["num_create"].sum() for p in range(4) for x in vc.loc[p]["num_create"]]
        vc["prop_alter_ref"] = [x / vc.loc[p]["num_alter_ref"].sum() for p in range(4) for x in vc.loc[p]["num_alter_ref"]]
        vc["prop_alter_alt"] = [x / vc.loc[p]["num_alter_alt"].sum() for p in range(4) for x in vc.loc[p]["num_alter_alt"]]
        vc["prop_destroy"] = [x / vc.loc[p]["num_destroy"].sum() for p in range(4) for x in vc.loc[p]["num_destroy"]]
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
    
    # save
    vc.to_csv(dst_vc, sep="\t")

    # Variable Amino Acids
    vc["phase"] = vc.index.get_level_values(0)
    vc.insert(loc=3, column="num_vc", value=1)
    vc.insert(loc=4, column="prop_vc", value=[1/len(vc.loc[p]) for p in range(4) for _ in range(len(vc.loc[p]))])
    vc_phase = [vc.loc[p] for p in range(4)]

    agg = {col : (lambda ser: ser.iloc[0]) for col in ["aattype"] + OUT}
    agg.update({col : "sum" for col in ["num_vc", "prop_vc"] + exp_col + freq_col})

    # DataFrame
    aat = pd.concat([v.groupby(["phase", "aa_ps", "aa_ds"]).aggregate(agg) for v in vc_phase])
    aat.to_csv(dst_aat, sep="\t")

    # Amino Acid Transition Types
    aat["phase"] = aat.index.get_level_values(0)
    aat.insert(loc=3, column="num_aat", value=1)
    aat.insert(loc=4, column="prop_aat", value=[1/len(aat.loc[p]) for p in range(4) for _ in range(len(aat.loc[p]))])
    aat_phase = [aat.loc[p] for p in [0,1,2]]

    agg = {col : "sum" for col in ["num_vc", "prop_vc", "num_aat", "prop_aat"] + exp_col + freq_col}

    # DataFrame    
    aatt = pd.concat([a.groupby(["phase", "aattype"]).aggregate(agg) for a in aat_phase])
    aatt.to_csv(dst_aatt, sep="\t")

# Count inserted/deleted amino acids
def count_outcomes(db:str=None):

    ''' Find the frequencies of amino acid outcomes
    
    **Source:** 
    * reference: `variants/proteome/nagnag_vc.txt`
    * variant: `variants/proteome/{db}_vc.txt`

    `variants/proteome/{nagnag/db}_outcomes.txt`
    --------------------------------------------
    * index: phase (0-3), cat ("ins", "del"), aa
    * columns:
        * reference: "exp_stoch", "exp_splice", "exp_as", "exp_ps", "exp_ds",
                     "num", "num_as", "num_ps", "num_ds", "prop", "prop_as", "prop_ps", "prop_ds"
        * variant: "exp_stoch", "exp_splice", 
                   "num", "num_create", "num_alter_ref", "num_alter_alt", "num_destroy",
                   "prop", "prop_create", "prop_alter_ref", "prop_alter_alt", "prop_destroy"
    '''

    # define {src}, {dst}, and {columns}
    if db:
        src = f"{ROOT}/variants/proteome/{db}_vc.txt"
        dst = f"{ROOT}/variants/proteome/{db}_outcomes.txt"

        columns = ["exp_stoch", "exp_splice", "num", "num_create", "num_alter_ref", "num_alter_alt", "num_destroy",
                   "prop", "prop_create", "prop_alter_ref", "prop_alter_alt", "prop_destroy"]
    else:
        src = f"{ROOT}/variants/proteome/nagnag_vc.txt"
        dst = f"{ROOT}/variants/proteome/nagnag_outcomes.txt"

        columns = ["exp_stoch", "exp_splice", "exp_as", "exp_ps", "exp_ds",
                   "num", "num_as", "num_ps", "num_ds", "prop", "prop_as", "prop_ps", "prop_ds"]

    # load variable codons
    vc = pd.read_csv(src, sep="\t")
    ph = vc["phase"]

    # inserted/deleted amino acids
    cats = list(itertools.product(OUT, AA))

    # count outcomes
    out = {p : {col : [vc[(ph == p) & (vc[c].str.contains(a, regex=False))][col].sum()
                       for c,a in cats] for col in columns} for p in range(4)}

    # DataFrame
    index = [(a,b,c) for a,(b,c) in itertools.product(range(4), cats)]
    df = pd.DataFrame({col : [x for p in range(4) for x in out[p][col]] for col in columns},
                      index=pd.MultiIndex.from_tuples(index, names=["phase", "cat", "aa"]))
    df.to_csv(dst, sep="\t")

# --- Combine Databases --- #

def consolidate_effects(suffix:str, n_index:int, n_cols:int, metric:str):

    ''' Combine reference and variant proteome effects
    
    **Source:** 
    * reference: `variants/proteome/nagnag_{suffix}.txt`
    * variant: `variants/proteome/{db}_{suffix}.txt`

    `variants/proteome/{metric}_freq.txt`
    -------------------------------------
    * rows: from `nagnag_{suffix}.txt` (first `n_index` columns)
    * columns: first `n_cols` columns of `nagnag_{suffix}.txt`,
    "exp_stoch", "exp_splice", "exp_ps", "exp_ds", "exp_as",
    "ref", "ref_ps", "ref_ds", "ref_as",
    "{db}_create", "{db}_alter_ref", "{db}_alter_alt", "{db}_destroy"

    `variants/proteome/{metric}_freq_count.txt`
    -------------------------------------
    * rows: from `nagnag_{suffix}.txt` (first `n_index` columns)
    * columns: first `n_cols` columns of `nagnag_{suffix}.txt`,
    "ref", "ref_ps", "ref_ds", "ref_as",
    "{db}_create", "{db}_alter_ref", "{db}_alter_alt", "{db}_destroy"
    '''

    # load frequencies
    ref = pd.read_csv(f"{ROOT}/variants/proteome/nagnag_{suffix}.txt" , sep="\t", index_col=list(range(n_index)))
    vnt = {db : pd.read_csv(f"{ROOT}/variants/proteome/{db}_{suffix}.txt", sep="\t", index_col=list(range(n_index))) for db in vnt_dbs}

    # DataFrames: count, proportions
    dfC = ref.iloc[:, :n_cols].copy()
    dfP = ref.iloc[:, :n_cols+5].copy()

    # reference frequencies
    for s in ["", "_ps", "_ds", "_as"]:
        dfC[f"ref{s}"] = ref[f"num{s}"]
        dfP[f"ref{s}"] = ref[f"prop{s}"]

    # variant frequencies
    for db,df in vnt.items():
        for s in ["", "_create", "_alter_ref", "_alter_alt", "_destroy"]:
            dfC[f"{db}{s}"] = df[f"num{s}"]
            dfP[f"{db}{s}"] = df[f"prop{s}"]

    # save
    dfP.to_csv(f"{ROOT}/variants/proteome/{metric}_freq.txt", sep="\t", index=True)
    dfC.to_csv(f"{ROOT}/variants/proteome/{metric}_freq_count.txt", sep="\t", index=True)

# --- Run --- #

def analyze_db(db:str=None):

    log_fn(db.replace("_", " ") if db else "hg38 reference", sub=True)

    if db:
        get_proteome_effects(db)
    
    count_transitions(db)
    count_outcomes(db)

# ===== RUN ===== #
log_script("11-analyze-variants.py")
log_fn("Calculating and counting variant proteomic effects")

analyze_db()
[analyze_db(db) for db in vnt_dbs]

# Consolidate
log_fn("Consolidating proteomic effects")
consolidate_effects(suffix="vc", n_index=3, n_cols=5, metric="variable_codons")
consolidate_effects(suffix="aat", n_index=3, n_cols=5, metric="variable_amino_acids")
consolidate_effects(suffix="aatt", n_index=2, n_cols=6, metric="amino_acid_transitions")
consolidate_effects(suffix="outcomes", n_index=3, n_cols=0, metric="amino_acid_outcomes")
