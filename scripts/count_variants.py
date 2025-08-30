''' Get frequencies:
* variant-affected splice sites (`variants/stats/vnt_site_freq.txt`)
* site-affecting variants (`variants/stats/vnt_freq.txt`)
* NAGNAG variant events + scenarios (`variants/stats/nagnag_vnt_event_scenario_freq.txt`)
* NAGNAG variant scenarios by phase (`variants/stats/nagnag_vnt_scenario_freq.txt`)
* variant base (`variants/stats/affected_base_freq.txt`)
* NAGNAG splice scenario phases (`stats/nagnag_scen_phase_freq.txt`)
'''

from utils import ROOT
from utils.lib import subprocess, progress_map, pd, np
from utils.seq import CHROMS
from utils.vnt import dtypes, vnt_dbs

def bed_wcl(filename:str) -> int:

    ''' Get length of .bed file (# of variants) '''

    output = subprocess.run(["wc", "-l", f"{ROOT}/variants/bed/{filename}.bed"], capture_output=True, text=True)
    return int(output.stdout.split()[0])

# Site frequencies
def count_sites():

    ''' Count variant-affected sites

    `variants/stats/vnt_site_freq.txt`
    ----------------------------------
    * columns: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice
    * rows: "3ss", "nagnag", "1off", "canon"
    '''

    print(" - variant-affected sites")

    # load sites
    ss3 = [pd.read_csv(f"{ROOT}/variants/affecting/{db}_3ss_vnt_containing_ssites.txt", sep="\t") for db in vnt_dbs]
    nn = [pd.read_csv(f"{ROOT}/variants/found/{db}_nagnag_vnt_containing_ssites.txt", sep="\t") for db in vnt_dbs]
    off = [ss[(ss["ssite_type"] == "1C") | (ss["ssite_type"] == "2NC")] for ss in nn]
    can = [ss[ss["ssite_type"] == "2C"] for ss in nn]

    # count possible
    n3 = pd.read_csv(f"{ROOT}/sites/3ss_uppercase.txt", sep="\t", dtype=dtypes).__len__()
    nO = pd.read_csv(f"{ROOT}/sites/1off_splice_sites.txt", sep="\t", dtype=dtypes).__len__()
    nC = pd.read_csv(f"{ROOT}/sites/nagnag_3ss.txt", sep="\t", dtype=dtypes).__len__()

    # DataFrames
    df = pd.DataFrame({
        "3ss" : [len(df) for df in ss3],
        "nagnag" : [len(df) for df in nn],
        "off" : [len(df) for df in off],
        "canon" : [len(df) for df in can]

    }, index=vnt_dbs).T
    df.insert(loc=0, column="total", value=[n3, nO+nC, nO, nC])

    df.to_csv(f"{ROOT}/variants/stats/vnt_site_freq.txt", sep="\t", index_label="database")

# Variant frequencies
def count_vnts():

    ''' Count site-affecting variants

    `variants/stats/vnt_freq.txt`
    -----------------------------
    * columns: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice
    * rows: "3ss", "nagnag", "create", "alter", "destroy"
    '''

    print(" - variants")

    # load variants
    ss3 = [pd.read_csv(f"{ROOT}/variants/affecting/{db}_3ss_vnts.txt", sep="\t", dtype=dtypes) for db in vnt_dbs]
    nn = [pd.read_csv(f"{ROOT}/variants/found/{db}_nagnag_affecting_vnts.txt", sep="\t", dtype=dtypes) for db in vnt_dbs]
    
    # count possible
    dbC = sum(progress_map(bed_wcl, [f"dbSNP_common_{c}" for c in CHROMS], n_cpu=24, disable=True))
    dbR = sum(progress_map(bed_wcl, [f"dbSNP_rare_{c}" for c in CHROMS], n_cpu=24, disable=True))
    cv = bed_wcl("ClinVar")
    hs = bed_wcl("HGMD_splice")

    # DataFrame
    df = pd.DataFrame({
        "total" : [dbC, dbR, cv, hs],
        "3ss" : [len(df) for df in ss3],
        "nagnag" : [len(df) for df in nn],
        "create" : [(df["vnt_effect"] == "CREATE").sum() for df in nn],
        "alter" : [(df["vnt_effect"] == "ALTER").sum() for df in nn],
        "destroy" : [(df["vnt_effect"] == "DESTROY").sum() for df in nn]
    }, index=vnt_dbs).T
    df.to_csv(f"{ROOT}/variants/stats/vnt_freq.txt", sep="\t", index_label="site")
    
# Event + scenario frequencies
def count_events_scens():

    ''' Count NAGNAG-affecting variant events and splice scenarios
    
    **Source:**
    * events: `vnt/found/{db}_nagnag_vnt_events.txt`
    * scenarios: `vnt/found/{db}_nagnag_vnt_scenarios.txt`

    `variants/stats/nagnag_vnt_event_scenario_freq.txt`
    ---------------------------------------------------
    * rows: (type, effect)
        * type: "events", "scenarios"
        * effect: "create", "alter", "destroy", "all"
    * columns: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"
    
    `variants/stats/nagnag_vnt_scenario_freq.txt`
    ---------------------------------------------
    * rows: (vnt_db, effect)
        * vnt_db: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"
        * effect: "create", "alter", "destroy", "all"
    * columns: "nc", "p0", "p1", "p2", "cds", "all"
    '''

    print(" - events and scenarios by variant effect")

    # load variant events + scenarios
    events = {db : pd.read_csv(f"{ROOT}/variants/found/{db}_nagnag_vnt_events.txt",
                               sep="\t", index_col=0, dtype=dtypes) for db in vnt_dbs}
    scens = {db : pd.read_csv(f"{ROOT}/variants/found/{db}_nagnag_vnt_scenarios.txt",
                              sep="\t", index_col=0, dtype=dtypes) for db in vnt_dbs}

    effects = ["CREATE", "ALTER", "DESTROY"]

    # events
    cE = np.array([[(df["vnt_effect"] == e).sum() for df in events.values()] for e in effects])
    cE = np.append(cE, np.sum(cE, axis=0).reshape(1, -1), axis=0)
    dfE = pd.DataFrame(cE, index=["create", "alter", "destroy", "all"], columns=vnt_dbs)

    # scenarios
    cS = np.array([[(df["vnt_effect"] == e).sum() for df in scens.values()] for e in effects])
    cS = np.append(cS, np.sum(cS, axis=0).reshape(1, -1), axis=0)
    dfS = pd.DataFrame(cS, index=["create", "alter", "destroy", "all"], columns=vnt_dbs)

    # combined
    df = pd.concat((dfE, dfS), axis=0, keys=["events", "scenarios"])
    df.to_csv(f"{ROOT}/variants/stats/nagnag_vnt_event_scenario_freq.txt",
              sep="\t", index_label=["type", "effect"])

    # by phase
    print(" - scenarios by phase")

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
        freq = pd.DataFrame(c, index=["create", "alter", "destroy", "all"],
                            columns=["nc", "p0", "p1", "p2", "cds", "all"])

        return freq

    # by phase
    dfs = [by_phase(db) for db in vnt_dbs]
    df = pd.concat(dfs, axis=0, keys=vnt_dbs)
    df.to_csv(f"{ROOT}/variants/stats/nagnag_vnt_scenario_freq.txt",
              sep="\t", index_label=["vnt_db", "effect"])

# Variant-affected base frequencies
def count_affected_bases():

    ''' Count NAGNAG-affecting variant bases
    
    src: variants/found/{db}_nagnag_vnt_scenarios.txt

    `variants/stats/affected_base_freq.txt`
    ----------=========--------------------
    rows: (effect, database)
    * effect: "CREATE", "ALTER", "DESTROY"
    * database: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice

    columns: "N1", "A1", "G1", "N2", "A2", "G2"
    '''

    print(" - variant-affected base")

    # load variant events
    events = [pd.read_csv(f"{ROOT}/variants/found/{db}_nagnag_vnt_events.txt", sep="\t", index_col=0, dtype=dtypes) for db in vnt_dbs]

    # possible positions for each variant effect type
    poss = {
        "CREATE" : [1, 2, 4, 5],
        "ALTER" : [0, 3],
        "DESTROY" : [1, 2, 4, 5]
    }

    # count frequencies
    counts = [[df[(df["vnt_effect"] == e) & (df["vnt_pos"] == p)].__len__()
               for df in events] for e,pos in poss.items() for p in pos]

    mpos = {0 : "N1", 1 : "A1", 2 : "G1", 3 : "N2", 4 : "A2", 5 : "G2"}
    [(e, mpos[p]) for e,pos in poss.items() for p in pos]


    # DataFrame
    df = pd.DataFrame(counts, index=pd.MultiIndex.from_tuples([(e, mpos[p])
                                                               for e,pos in poss.items() for p in pos],
                                                               names=["effect", "base"]),
                      columns=vnt_dbs)
    df.to_csv(f"{ROOT}/variants/stats/affected_base_freq.txt", sep="\t")

print("computing variant frequency statistics...")
count_sites()
count_vnts()
count_events_scens()
count_affected_bases()
print("done")

# dbSNP Common: 33,629,539
# dbSNP Rare: 56,354,5478
# ClinVar: 3401,768
# HGMD: 35,462
