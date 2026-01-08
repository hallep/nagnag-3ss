''' Get frequencies:
* variant-affected splice sites (`variants/stats/vnt_site_freq.txt`)
* site-affecting variants (`variants/stats/vnt_freq.txt`)
* NAGNAG variant events + scenarios (`variants/stats/nagnag_vnt_event_scenario_freq.txt`)
* NAGNAG variant scenarios by phase (`variants/stats/nagnag_vnt_scenario_freq.txt`)
* variant base (`variants/stats/affected_base_freq.txt`)
* NAGNAG splice scenario phases (`stats/nagnag_scen_phase_freq.txt`)
'''

from lib import argparse, subprocess, pd, np
from utils import ROOT, log_script, log_fn, single_map
from utils.seq import CHROMS
from utils.vnt import dtypes, all_vnt_dbs

# Variant Databases
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--num-threads", type=int, default=24, help="Number of parallel threads available")
parser.add_argument("-H", "--ignore-HGMD", action="store_true", help="Do not process/analyze HGMD Splice variants")
args = parser.parse_args()

vnt_dbs = all_vnt_dbs
if args.ignore_HGMD:
    vnt_dbs.remove("HGMD_splice")

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

    }, index=pd.Index(vnt_dbs, names=["db"])).T
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

    # load variants
    ss3 = [pd.read_csv(f"{ROOT}/variants/affecting/{db}_3ss_vnts.txt", sep="\t", dtype=dtypes) for db in vnt_dbs]
    nn = [pd.read_csv(f"{ROOT}/variants/found/{db}_nagnag_affecting_vnts.txt", sep="\t", dtype=dtypes) for db in vnt_dbs]
    
    # count possible
    dbC = sum(single_map(bed_wcl, [f"dbSNP_common_{c}" for c in CHROMS], n_procs=args.num_threads))
    dbR = sum(single_map(bed_wcl, [f"dbSNP_rare_{c}" for c in CHROMS], n_procs=args.num_threads))
    cv = bed_wcl("ClinVar")

    if not args.ignore_HGMD:
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
    
    `variants/stats/nagnag_vnt_scenario_phase_freq.txt`
    ---------------------------------------------
    * rows: (vnt_db, effect)
        * vnt_db: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"
        * effect: "create", "alter", "destroy", "all"
    * columns: "nc", "p0", "p1", "p2", "cds", "all"

    `variants/stats/nagnag_vnt_scenario_stype_freq.txt`
    ---------------------------------------------
    * rows: (vnt_db, effect)
        * vnt_db: "dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"
        * effect: "create", "alter", "destroy", "all"
    * columns: "ps", "ds", "as", "all"
    '''

    # load variant events + scenarios
    events = {db : pd.read_csv(f"{ROOT}/variants/found/{db}_nagnag_vnt_events.txt",
                               sep="\t", index_col=0, dtype=dtypes) for db in vnt_dbs}
    scens = {db : pd.read_csv(f"{ROOT}/variants/found/{db}_nagnag_vnt_scenarios.txt",
                              sep="\t", index_col=0, dtype=dtypes) for db in vnt_dbs}

    effects = ["CREATE", "ALTER", "DESTROY"]

    # events
    cE = np.array([[(df["vnt_effect"] == e).sum() for df in events.values()] for e in effects])
    cE = np.append(cE, np.sum(cE, axis=0).reshape(1, -1), axis=0)
    dfE = pd.DataFrame(cE, index=pd.Index(["create", "alter", "destroy", "all"], name="effect"), columns=vnt_dbs)

    # scenarios
    cS = np.array([[(df["vnt_effect"] == e).sum() for df in scens.values()] for e in effects])
    cS = np.append(cS, np.sum(cS, axis=0).reshape(1, -1), axis=0)
    dfS = pd.DataFrame(cS, index=pd.Index(["create", "alter", "destroy", "all"], name="effect"), columns=vnt_dbs)

    # combined
    df = pd.concat((dfE, dfS), axis=0, keys=["events", "scenarios"])
    df.to_csv(f"{ROOT}/variants/stats/nagnag_vnt_event_scenario_freq.txt",
              sep="\t", index_label=["type", "effect"])

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
        freq = pd.DataFrame(c, index=pd.Index(["create", "alter", "destroy", "all"], name="effect"),
                            columns=["nc", "p0", "p1", "p2", "cds", "all"])

        return freq

    # by phase
    dfs = [by_phase(db) for db in vnt_dbs]
    df = pd.concat(dfs, axis=0, keys=vnt_dbs)
    df.to_csv(f"{ROOT}/variants/stats/nagnag_vnt_scenario_phase_freq.txt",
              sep="\t", index_label=["vnt_db", "effect"])

    # by scenario splice type
    effects = ["CREATE", "ALTER", "DESTROY"]
    stypes = ["PS", "DS", "AS"]

    def by_stype(db:str) -> pd.DataFrame:

        df = scens[db]

        # count by phase + effect
        c = [[df[(df["splice_type"] == s) & (df["vnt_effect"] == e)].__len__() for s in stypes] for e in effects]
        c = np.append(c, np.sum(c, axis=0).reshape(1, -1), axis=0)
        c = np.append(c, np.sum(c, axis=1).reshape(-1, 1), axis=1)

        # create DataFrame
        freq = pd.DataFrame(c, index=pd.Index(["create", "alter", "destroy", "all"], name="effect"),
                            columns=["ps", "ds", "as", "all"])

        return freq

    # by phase
    dfs = [by_stype(db) for db in vnt_dbs]
    df = pd.concat(dfs, axis=0, keys=vnt_dbs, names=["db"])
    df.to_csv(f"{ROOT}/variants/stats/nagnag_vnt_scenario_stype_freq.txt", sep="\t", index_label=["vnt_db", "effect"])

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

    # load variant events
    events = {db : pd.read_csv(f"{ROOT}/variants/found/{db}_nagnag_vnt_events.txt", sep="\t", index_col=0, dtype=dtypes) for db in vnt_dbs}

    # possible positions for each variant effect type
    poss = {
        "CREATE" : [1, 2, 4, 5],
        "ALTER" : [0, 3],
        "DESTROY" : [1, 2, 4, 5]
    }
    mpos = {0 : "N1", 1 : "A1", 2 : "G1", 3 : "N2", 4 : "A2", 5 : "G2"}

    # count frequencies for each database by effect
    count = {
        e : pd.DataFrame({
            db : [df[(df["vnt_effect"] == e) & (df["vnt_pos"] == p)].__len__() for p in pos] for db,df in events.items()
        }, index=pd.Index([mpos[p] for p in pos], name="base")) for e,pos in poss.items()
    }

    # calculate proportions
    prop = {e : df / df.sum() for e,df in count.items()}

    # combine counts and proportions
    freq = {e : pd.concat((count[e], prop[e]), axis=1, keys=["num", "prop"]) for e in list(poss)}

    # combine effects
    df = pd.concat((freq.values()), keys=freq.keys(), names=["effect"])
    df.to_csv(f"{ROOT}/variants/stats/affected_base_freq.txt", sep="\t")

# ===== RUN ===== #
if __name__ == "__main__":
    # log_script("10-count-variants.py")
    # log_fn("Computing variant frequency statistics")

    # log_fn("variant-affected sites", sub=1)
    # count_sites()

    # log_fn("variants", sub=1)
    # count_vnts()

    log_fn("events and scenarios by variant effect", sub=1)
    count_events_scens()

    # log_fn("variant-affected base", sub=1)
    # count_affected_bases()
