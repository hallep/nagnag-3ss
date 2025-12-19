''' Get site frequencies:
* splice site types (`stats/3ss_site_type_freq.txt`)
* NAGNAGs, expected vs. observed (`stats/exp_obs_nagnag_freq.txt`)
* 1-NAG motifs (`stats/1nag_motif_freq.txt`)
* NAGNAG motifs (`stats/nagnag_motif_freq.txt`)
* NAGNAG splice types (`stats/nagnag_splice_type_freq.txt`)
* NAGNAG splice scenario phases (`stats/nagnag_scen_phase_freq.txt`)
'''

from lib import itertools, pd, np, proportions_ztest
from utils import ROOT, log_script, log_fn
from utils.seq import N

# Splice site types
def ss_type_freq():

    ''' Determine the frequency of each splice site type (amongst uppercase splice sites)
    
    **Source:** `sites/3ss.txt`

    `stats/3ss_site_type_freq.txt`
    ------------------------------
    * **sstype** (*str*): splice site type
        * possible values: "1C", "1NC", "2C", "2NC", "3+", "all"
    * **num_ss** (*int*): number of splice sites of each type
    * **prop_ss** (*float*): proportion of splice sites of each type
    * **num_cs** (*int*): number of cleavage sites of each type
    * **prop_cs** (*float*): proportion of cleavage sites of each type
    * **num_scen** (*int*): int number of splice scenarios of each type
    * **prop_scen** (*float*): proportion of splice scenarios of each type
    '''
    
    # load splice sites
    ssites = pd.read_csv(f"{ROOT}/sites/3ss.txt", sep="\t", index_col=0)
    ssites = ssites[ssites["uppercase"] == 1]

    sstypes = ["1C", "1NC", "2C", "2NC", "3+"]

    # count splice sites
    count_ss = np.array([(ssites["ssite_type"] == s).sum() for s in sstypes])
    count_ss = np.append(count_ss, sum(count_ss))

    # count cleavage sites
    count_cs = np.array([ssites[ssites["ssite_type"] == s]["num_csites"].sum() for s in sstypes])
    count_cs = np.append(count_cs, sum(count_cs))

    # count splice scenarios
    count_scen = np.array([ssites[ssites["ssite_type"] == s]["num_scens"].sum() for s in sstypes])
    count_scen = np.append(count_scen, sum(count_scen))

    # DataFrame
    df = pd.DataFrame({
        "sstype" : sstypes + ["all"],
        "num_ss" : count_ss,
        "prop_ss" : count_ss / sum(count_ss[:-1]),
        "num_cs" : count_cs,
        "prop_cs" : count_cs / sum(count_cs[:-1]),
        "num_scen" : count_scen,
        "prop_scen" : count_scen / sum(count_scen[:-1])
    }).set_index(keys="sstype", inplace=False)
    df.to_csv(f"{ROOT}/stats/3ss_site_type_freq.txt", sep="\t")

# Expected vs. observed NAGNAG frequency
def exp_obs_nagnag_freq():

    ''' Compare the stochastic expected and observed frequency of NAGNAGs
    
    **Source:** `sites/3ss.txt`
    
    `stats/exp_obs_nagnag_freq.txt`
    -------------------------------
    * **exp** (*float*): expected frequency of NAGNAGs under stochastic model
    * **obs** (*float*): proportion of NAG-containing 3' splice sites that are NAGNAGs
    * **x** (*int*): number of NAGNAG 3' splice sites
    * **n** (*int*): number of NAG-containing 3' splice sites
    * **z** (*float*): z-score from 1-proportion 2-sided z-test
    * **p** (*float*): p-value in 1-proportion 2-sided z-test
    '''

    # load splice sites
    ssites = pd.read_csv(f"{ROOT}/sites/3ss.txt", sep="\t", index_col=0)
    ssites = ssites[(ssites["uppercase"] == 1)]

    num_nagnag = int((ssites["ssite_type"] == "2C").sum())
    num_nag = int(ssites["has_nag"].sum())
    obs = num_nagnag / num_nag

    n = (1/4) * (1/4)
    nn = 1 - n
    exp = 2 * (nn * n * nn)

    # z-test
    z,p = proportions_ztest(count=num_nagnag, nobs=num_nag, value=exp)

    # Series
    ser = pd.Series(data=[exp, obs, num_nagnag, num_nag, z, p],
                    index=["exp", "obs", "x", "n", "z", "p"], name="exp_obs")
    ser.to_csv(f"{ROOT}/stats/exp_obs_nagnag_freq.txt", sep="\t")

# 1-NAG motifs
def canon_1nag_motif_freq():

    ''' Determine the frequency of (uppercase) canonical 1-NAG motifs
    
    **Source:** `sites/3ss.txt`
    
    `stats/1nag_motif_freq.txt`
    ---------------------------
    * **motif** (*str*): 1-NAG motif
        * possible values: "AAG", "CAG", "GAG", "TAG", "all"
    * **num_ss** (*int*): number of 1-NAG splice sites of each motif
    * **prop_ss** (*float*): proportion of 1-NAG splice sites of each motif
    * **num_scen** (*int*): number of 1-NAG splice scenarios of each motif
    * **prop_scen** (*float*): proportion of 1-NAG splice scenarios of each motif
    '''

    # load splice sites
    ssites = pd.read_csv(f"{ROOT}/sites/3ss.txt", sep="\t", index_col=0)
    ssites = ssites[(ssites["uppercase"] == 1) & (ssites["ssite_type"] == "1C")]

    motifs = ["AAG", "CAG", "GAG", "TAG"]

    # count splice sites
    count_ss = np.array([(ssites["ssite_seq"] == m).sum() for m in motifs])
    count_ss = np.append(count_ss, sum(count_ss))

    # count splice scenarios
    count_scen = np.array([ssites[ssites["ssite_seq"] == m]["num_scens"].sum() for m in motifs])
    count_scen = np.append(count_scen, sum(count_scen))

    # DataFrame
    df = pd.DataFrame({
        "motif" : motifs + ["all"],
        "num_ss" : count_ss,
        "prop_ss" : count_ss / sum(count_ss[:-1]),
        "num_scen" : count_scen,
        "prop_scen" : count_scen / sum(count_scen[:-1])
    }).set_index(keys="motif", inplace=False)
    df.to_csv(f"{ROOT}/stats/1nag_motif_freq.txt", sep="\t")

# NAGNAG motifs
def nagnag_motif_freq():

    ''' Determine the frequency of NAGNAG motifs
    
    **Source:** `sites/nagnag_3ss.txt`
    
    `stats/nagnag_motif_freq.txt`
    ---------------------------
    * **motif** (*str*): NAGNAG motif
    * **num_all**, **num_ps**, **num_ds**, **num_as** (*int*): 
        number of all, proximally-, distally-, and alternatively-spliced NAGNAGs (respectively) of each motif
    * **prop_all**, **prop_ps**, **prop_ds**, **prop_as** (*float*): 
        proportion of all, proximally-, distally-, and alternatively-spliced NAGNAGs (respectively) of each motif
    '''

    # load splice sites
    ssites = pd.read_csv(f"{ROOT}/sites/nagnag_3ss.txt", sep="\t", index_col=0, dtype={"csite_pos":"str"})
    
    cs_pos = [ssites["csite_pos"], "0", "1", "0,1"]
    motifs = [f"{n1}AG{n2}AG" for n1,n2 in itertools.product(N, repeat=2)]

    # count splice sites
    freqs = np.array([[len(ssites[(ssites["csite_pos"] == p) & (ssites["ssite_seq"] == m)]) for m in motifs] for p in cs_pos])
    freqs = np.concatenate((freqs, np.sum(freqs, axis=1, keepdims=True)), axis=1)
    
    suffix = ["all", "ps", "ds", "as"]
    d = {"motif" : motifs + ["all"]}
    d.update({f"num_{s}" : f for s,f in zip(suffix, freqs)})
    d.update({f"prop_{s}" : f / sum(f[:-1]) for s,f in zip(suffix, freqs)})

    # DataFrame
    df = pd.DataFrame(d).set_index(keys="motif", inplace=False)
    df.to_csv(f"{ROOT}/stats/nagnag_motif_freq.txt", sep="\t")

# NAGNAG splice types
def nagnag_splice_type_freq():
    
    ''' Determine the frequency of NAGNAG splice types
    
    **Source:**
    * `sites/nagnag_3ss.txt`
    * `sites/nagnag_scens.txt`
    
    `stats/nagnag_splice_type_freq.txt`
    -----------------------------------
    * **stype** (*str*): NAGNAG splice type
        * possible values: "PS", "DS", "AS", "all"
    * **num_ss** (*int*): number of NAGNAG splice sites of each splice type
    * **prop_ss** (*float*): proportion of NAGNAG splice sites of each splice type
    * **num_scen** (*int*): number of NAGNAG splice scenarios of each splice type
    * **prop_scen** (*float*): proportion of NAGNAG splice scenarios of each splice type
    '''

    # load splice sites and scenarios
    ssites = pd.read_csv(f"{ROOT}/sites/nagnag_3ss.txt", sep="\t", index_col=0)
    scens = pd.read_csv(f"{ROOT}/sites/nagnag_scens.txt", sep="\t", index_col=0)

    pos = ["0", "1", "0,1"]
    stypes = ["PS", "DS", "AS"]

    # count splice sites
    count_ss = np.array([(ssites["csite_pos"] == p).sum() for p in pos])
    count_ss = np.append(count_ss, sum(count_ss))

    # count splice scenarios
    count_scen = np.array([(scens["splice_type"] == s).sum() for s in stypes])
    count_scen = np.append(count_scen, sum(count_scen))

    # DataFrame
    df = pd.DataFrame({
        "stype" : ["PS", "DS", "AS", "all"],
        "num_ss" : count_ss,
        "prop_ss" : count_ss / sum(count_ss[:-1]),
        "num_scen" : count_scen,
        "prop_scen" : count_scen / sum(count_scen[:-1])
    }).set_index(keys="stype", inplace=False)
    df.to_csv(f"{ROOT}/stats/nagnag_splice_type_freq.txt", sep="\t")

# NAGNAG splice scenario phases
def nagnag_scen_phase_freq():
    
    ''' Determine the frequency of NAGNAG splice scenarios by phase
    
    **Source:** `sites/nagnag_scens.txt`
    
    `stats/nagnag_scen_phase_freq.txt`
    ----------------------------------
    * **phase** (*int*): phase of downstream exon
        * possible values: -1, 0, 1, 2, 3
        * -1 indicates non-coding
        * 3 indicates total
    * **num** (*int*): number of NAGNAG splice scenarios in each phase
    * **prop** (*float*): proportion of NAGNAG splice scenarios in each phase
    * **prop_CDS** (*float*): proportion of coding sequence NAGNAG splice scenarios in each phase
    '''

    # load splice scenarios
    scens = pd.read_csv(f"{ROOT}/sites/nagnag_scens.txt", sep="\t", index_col=0)

    phase = [-1, 0, 1, 2]
    
    # get counts
    count = np.array([(scens["phase"] == p).sum() for p in phase])
    count = np.append(count, sum(count))

    # DataFrame
    df = pd.DataFrame({
        "phase" : [-1, 0, 1, 2, 3],
        "num" : count,
        "prop" : count / sum(count[:-1]),
        "prop_CDS" : [np.NaN] + list(count[1:-1] / sum(count[1:-1])) + [1]
    }).set_index(keys="phase", inplace=False)
    df.to_csv("stats/nagnag_scen_phase_freq.txt", sep="\t")

# ===== RUN ===== #
log_script("01-sites.py")
log_fn("Computing splice site frequency statistics")

log_fn("3' splice site types", sub=True)
ss_type_freq()
log_fn("Expected NAGNAG frequency", sub=True)
exp_obs_nagnag_freq()

# Motifs
log_fn("1-NAG motifs", sub=True)
canon_1nag_motif_freq()
log_fn("NAGNAG motifs", sub=True)
nagnag_motif_freq()

log_fn("NAGNAG splice types", sub=True)
nagnag_splice_type_freq()
log_fn("NAGNAG splice scenario phases", sub=True)
nagnag_scen_phase_freq()
