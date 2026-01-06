''' Get position-wise base frequencies for 1-NAG and NAGNAG 3' splice sites '''

from lib import pd, np
from utils import ROOT, log_script, log_fn
from utils.seq import N

# 1-NAG position-wise base frequencies
def poswise_freq_1nag():

    ''' Calculate position-wise frequencies of upstream and downstream exon bases of 1-NAG splice scenarios

    **Source:** `sites/3cs.txt`

    `proteome/1nag_poswise_freq.txt`
    -----------------------------
    **syntax**: df[phase][splice type][position][base]
    * **phase**: "all", "nc", "p0", "p1", "p2", "cds"
    * **splice type**: "all"
    * **position**: "u-2", "u-1", "d1", "d2"
    * **base**: "A", "C", "G", "T"

    `proteome/1nag_poswise_freq_count.txt`
    ----------------------------------------
    same as `proteome/1nag_poswise_freq.txt`, but with raw counts
    '''

    # load 1-NAG cleavage sites
    n1 = pd.read_csv(f"{ROOT}/sites/3cs.txt", sep="\t", index_col=0)
    n1 = n1[(n1["ssite_type"] == "1C") & (n1["csite_seq"].str.isupper())]

    # split into splice scenarios
    phases = n1["phase"]
    up_seqs = n1["scen_eflank_5ss"]
    down_seqs = n1["csite_eflank_3ss"]
    num_scens = n1["num_scens"]

    scens = pd.DataFrame({
        "u-2" : [u[-2].upper() for up in up_seqs for u in up.split(",")],
        "u-1" : [u[-1].upper() for up in up_seqs for u in up.split(",")],
        "d1" : [d[0].upper() for n,d in zip(num_scens, down_seqs) for _ in range(n)],
        "d2" : [d[1].upper() for n,d in zip(num_scens, down_seqs) for _ in range(n)],
        "p" : [int(p) for phase in phases for p in phase.split(",")]
    })

    # position-wise frequencies
    phase = [scens["p"], -1, 0, 1, 2]
    position = ["u-2", "u-1", "d1", "d2"]

    freq = np.array([[[len(scens[(scens[pos] == n) & (scens["p"] == p)]) for n in N] for pos in position] for p in phase])
    freq = np.concatenate((freq, np.sum(freq[-3:], axis=0, keepdims=True)))

    # DataFrames
    phase = ["all", "nc", "p0", "p1", "p2", "cds"]

    # counts
    count = pd.DataFrame({(ph, "all", pos) : freq[i][j] for i,ph in enumerate(phase)
                          for j,pos in enumerate(position)}, index=pd.Index(N, name="base"))
    count.to_csv(f"{ROOT}/proteome/1nag_poswise_freq_count.txt", sep="\t", index=True)

    # proportions
    prop = pd.DataFrame({(ph, "all", pos) : freq[i][j]/sum(freq[i][j]) for i,ph in enumerate(phase)
                         for j,pos in enumerate(position)}, index=pd.Index(N, name="base"))
    prop.to_csv(f"{ROOT}/proteome/1nag_poswise_freq.txt", sep="\t", index=True)

# NAGNAG position-wise base frequencies
def poswise_freq_nagnag():

    ''' Calculate position-wise frequencies of upstream and downstream exon bases of NAGNAG splice scenarios

    **Source:** `sites/nagnag_scens.txt`

    `proteome/nagnag_poswise_freq.txt`
    ----------------------------------
    **syntax**: df[phase][splice type][position][base]
    * **phase**: "all", "nc", "p0", "p1", "p2", "cds"
    * **splice type**: "all", AS", "PS", "DS"
    * **position**: "u-2", "u-1", "d1", "d2"
    * **base**: "A", "C", "G", "T"

    `proteome/nagnag_poswise_freq_count.txt`
    ----------------------------------------
    same as `proteome/nagnag_poswise_freq.txt`, but with raw counts
    '''
    
    # load splice scenarios
    n2 = pd.read_csv(f"{ROOT}/sites/nagnag_scens.txt", sep="\t", index_col=0)

    scens = pd.DataFrame({
        "u-2" : n2["up_seq"].str.upper().str.get(-2),
        "u-1" : n2["up_seq"].str.upper().str.get(-1),
        "n1" : n2["ssite_seq"].str.upper().str.get(0),
        "n2" : n2["ssite_seq"].str.upper().str.get(3),
        "d1" : n2["down_seq"].str.upper().str.get(0),
        "d2" : n2["down_seq"].str.upper().str.get(1),
        "p" : n2["phase"],
        "s" : n2["splice_type"]
    })

    phase = [scens["p"], -1, 0, 1, 2]
    stype = [scens["s"], "AS", "PS", "DS"]
    position = ["u-2", "u-1", "n1", "n2", "d1", "d2"]

    # position-wise frequencies
    freq = np.array([[[[len(scens[(scens["p"] == ph) & (scens["s"] == st) & (scens[pos] == n)])
                        for n in N] for pos in position] for st in stype] for ph in phase])
    freq = np.concatenate((freq, np.sum(freq[-3:], axis=0, keepdims=True)))
    
    # DataFrame
    phase = ["all", "nc", "p0", "p1", "p2", "cds"]
    stype = ["all", "AS", "PS", "DS"]

    # counts
    count = pd.DataFrame({(ph, st, pos) : freq[i][j][k] for i,ph in enumerate(phase)
                          for j,st in enumerate(stype) for k,pos in enumerate(position)}, index=pd.Index(N, name="base"))
    count.to_csv(f"{ROOT}/proteome/nagnag_poswise_freq_count.txt", sep="\t", index=True)

    # proportions
    prop = pd.DataFrame({(ph, st, pos) : freq[i][j][k]/sum(freq[i][j][k]) for i,ph in enumerate(phase)
                         for j,st in enumerate(stype) for k,pos in enumerate(position)}, index=pd.Index(N, name="base"))
    prop.to_csv(f"{ROOT}/proteome/nagnag_poswise_freq.txt", sep="\t", index=True)

# ===== RUN ===== #
log_script("06-poswise.py")
log_fn("calculating position-wise base frequencies")

log_fn("1-NAGs", sub=1)
poswise_freq_1nag()

log_fn("NAGNAGs", sub=1)
poswise_freq_nagnag()
