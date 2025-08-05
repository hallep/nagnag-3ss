''' Get all possible NAGNAG variable codons

Source:
-------
* **splice scenarios**: `sites/nagnag_scens.txt`
* **1-NAGs**: `proteome/1nag_poswise_freq.txt`
* **NAGNAGs**: `proteome/nagnag_poswise_freq.txt`

`proteome/variable_codon_freq.txt`
-------------------------------
* index: (phase, vc_ps, vc_ds)
* columns:
    * ("tsn", x); x = {"aa_ps", "aa_ds", "aattype", "ins_aa", "del_aa"}
    * ("exp", x); x = {"exp_stoch", "exp_splice", "exp_ps", "exp_ds", "exp_as"}
    * (freq, stype); freq = {"prop", "num"}; stype = {"all", "ps", "ds", "as"}

`proteome/variable_aa_freq.txt`
----------------------------
* index: (phase, aa_ps, aa_ds)
* columns:
    * ("tsn", x); x = {"aattype", "ins_aa", "del_aa"}
    * ("exp", x); x = {"exp_stoch", "exp_splice", "exp_ps", "exp_ds", "exp_as"}
    * (freq, stype); freq = {"prop", "num"}; stype = {"all", "ps", "ds", "as"}

`proteome/transition_type.txt`
---------------------------
* index: (phase, aattype)
* columns:
    * ("exp", x); x = {"exp_stoch", "exp_splice", "exp_ps", "exp_ds", "exp_as"}
    * (freq, stype); freq = {"prop", "num"}; stype = {"all", "ps", "ds", "as"}
'''

from lib import ROOT, itertools, pd, np
from sequence import N, AA, get_tsn, get_categorical_outcomes

# Expected variable codon frequency
def exp_vc_freq(phase:int, wsource:int=0, stype:str="all") -> list[float]:

    ''' Get expected variable codon frequencies

    Parameters
    ----------
    phase : int {0, 1, 2}
        NAGNAG phase
    wsource : int {0, 1}
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

    # position-wise base frequencies
    f = {
        # stochastic
        0 : pd.DataFrame(np.full((4, (5*4*6)), 1/4), index=N,
                         columns=pd.MultiIndex.from_product((["all", "nc", "p0", "p1", "p2"],
                                                             ["all", "AS", "PS", "DS"],
                                                             ["u-2", "u-1", "n1", "n2", "d1", "d2"]))),

        # weighted from 1-NAGs
        1 : pd.read_csv(f"{ROOT}/proteome/1nag_poswise_freq.txt", sep="\t", index_col=0, header=[0,1,2]),

        # weighted from NAGNAGs        
        2 : pd.read_csv(f"{ROOT}/proteome/nagnag_poswise_freq.txt", sep="\t", index_col=0, header=[0,1,2])
    }

    # base positions
    pos = {
        0 : ["n2"],
        1 : ["u-1", "n2", "d1", "d2"],
        2 : ["u-2", "u-1", "n2", "d1"]
    }[phase]

    # weight sources
    wgts = [2 if (x == "n2") else 1 for x in pos] if (wsource == 1) else [0] * len(pos)

    # phases
    phases = ["nc" if (x == "n2") else f"p{phase}" for x in pos]

    # splice types
    stypes = [stype if (x == "n2") else "all" for x in pos]

    # variable bases
    r = 1 if (phase == 0) else 4
    variable = list(itertools.product(N, repeat=r))

    return [np.prod([f[w][p][s][x][b] for w,p,s,x,b in zip(wgts, phases, stypes, pos, vb)]) for vb in variable]

# load scenarios
scens = pd.read_csv(f"{ROOT}/sites/nagnag_scens.txt", sep="\t")
scens = scens[scens["phase"] != -1]
ph = scens["phase"]
st = scens["splice_type"]
ps = scens["vc_ps"].str.upper()
ds = scens["vc_ds"].str.upper()

stypes = [st, "PS", "DS", "AS"]
st_cols = ["all", "ps", "ds", "as"]

# Variable Codons
def get_vc(phase:int) -> pd.DataFrame:
    
    # variable bases
    vb = list(itertools.product(N, repeat=4))

    # amino acid transition
    if phase == 0:
        tsn = [get_tsn(up="NNN", motif=f"NAG{n}AG", down="NNN", phase=0) for n in N]
    elif phase == 1:
        tsn = [get_tsn(up=f"NN{v[0]}", motif=f"NNN{v[1]}AG", down=f"{v[2]}{v[3]}N", phase=1) for v in vb]
    elif phase == 2:
        tsn = [get_tsn(up=f"N{v[0]}{v[1]}", motif=f"NNN{v[2]}AG", down=f"{v[3]}NN", phase=2) for v in vb]

    vc_ps, vc_ds, aa_ps, aa_ds, aattype = zip(*tsn)
    ins_aa, del_aa = zip(*[get_categorical_outcomes(1, p, d) for p,d in zip(aa_ps, aa_ds)])
    freq = np.array([[len(scens[(ph == phase) & (st == s) & (ps == p) & (ds == d)]) for p,d in zip(vc_ps, vc_ds)] for s in stypes])

    d = {
        ("tsn", "aa_ps") : aa_ps,
        ("tsn", "aa_ds") : aa_ds,
        ("tsn", "aattype") : aattype,
        ("tsn", "ins_aa") : ins_aa,
        ("tsn", "del_aa") : del_aa,
        ("exp", "exp_stoch") : exp_vc_freq(phase=phase, wsource=0, stype="all"),
        ("exp", "exp_splice") : exp_vc_freq(phase=phase, wsource=1, stype="all"),
        ("exp", "exp_ps") : exp_vc_freq(phase=phase, wsource=1, stype="PS"),
        ("exp", "exp_ds") : exp_vc_freq(phase=phase, wsource=1, stype="DS"),
        ("exp", "exp_as") : exp_vc_freq(phase=phase, wsource=1, stype="AS"),
    }
    d.update({("num", s) : f for s,f in zip(st_cols, freq)})
    d.update({("prop", s) : f/sum(f) for s,f in zip(st_cols, freq)})

    vc = pd.DataFrame(d).sort_values(by=("exp", "exp_splice"), ascending=False, inplace=False)
    vc.index = pd.MultiIndex.from_tuples(zip(vc_ps, vc_ds))

    return vc

# by phase
vcs = {p : get_vc(p) for p in range(3)}
vc = pd.concat(vcs.values(), keys=vcs.keys())
vc.index.names = ["phase", "vc_ps", "vc_ds"]
vc.to_csv(f"{ROOT}/proteome/variable_codon_freq.txt", sep="\t")

# Variable Amino Acids
def get_aas(vc:pd.DataFrame) -> pd.DataFrame:
    agg_fn = {c : "sum" if (c[0] != "tsn") else lambda x: x.iloc[0] for c in vc.columns.values[2:]}
    aa = vc.groupby([("tsn", "aa_ps"), ("tsn", "aa_ds")]).aggregate(agg_fn)
    aa.index.names = ["aa_ps", "aa_ds"]
    aa.sort_values(by=("exp", "exp_splice"), ascending=False, inplace=True)
    return aa

# by phase
aas = {p : get_aas(a) for p,a in vcs.items()}
aa = pd.concat(aas.values(), keys=aas.keys())
aa.index.names = ["phase", "aa_ps", "aa_ds"]
aa.to_csv(f"{ROOT}/proteome/variable_aa_freq.txt", sep="\t")

# Amino Acid Transition Type
def get_aatt(aa:pd.DataFrame) -> pd.DataFrame:
    agg_fn = {c : "sum" for c in aa.columns.values[3:]}
    aattype = aa.groupby(("tsn", "aattype")).aggregate(agg_fn)
    aattype.index.names = ["aattype"]
    aattype.sort_values(by=("exp", "exp_splice"), ascending=False, inplace=True)
    return aattype

# by phase
aatts = {p : get_aatt(a) for p,a in aas.items()}
aatt = pd.concat(aatts.values(), keys=aatts.keys())
aatt.index.names = ["phase", "aattype"]
aatt.to_csv(f"{ROOT}/proteome/transition_type.txt", sep="\t")

# All Scenarios
def combine(dfs:dict[int, pd.DataFrame]) -> pd.DataFrame:
    df = pd.concat(dfs.values())
    df.index.names = dfs[0].index.names

    # expected
    freqs = {p : (scens["phase"] == p).sum() / len(scens) for p in range(3)}    
    for col in itertools.product(["exp"], ["exp_stoch", "exp_splice", "exp_ps", "exp_ds", "exp_as"]):        
        df[col] = pd.concat([dfs[p][col] * freqs[p] for p in range(3)]).values

    # observed proportions
    for c in ["all", "ps", "ds", "as"]:
        df[("prop", c)] = df[("num", c)] / len(scens)
    
    return df

vcs[3] = combine(vcs)

# Amino Acid Outcomes
def get_out(vc:pd.DataFrame) -> pd.DataFrame:

    def cat(col:str) -> pd.DataFrame:
        return pd.DataFrame({c : [vc[vc[("tsn", col)].str.contains(a, regex=False)][c].sum() for a in AA]
                                for c in vc.columns.values[6:]}, index=AA)
    
    out = pd.concat([cat("ins_aa"), cat("del_aa")], keys=["ins", "del"])
    out.index.names = ["cat", "aa"]
    return out

outs = {p : get_out(v) for p,v in vcs.items()}
out = pd.concat(outs.values(), keys=outs.keys())
out.index.names = ["phase", "cat", "aa"]
out.to_csv(f"{ROOT}/proteome/aa_outcomes.txt", sep="\t")
