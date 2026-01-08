''' Assess the variant effect 

`variants/found/{db}_nagnag_vnt_events.txt`
------------------------
* **dStrength** (*str*): change in NAG strength due to variant
    * NAGNAG-creating: [created] - [original]
    * NAGNAG-altering: [alternate] - [reference]
    * NAGNAG-destroying: [destroyed] - [unaffected]
'''

from lib import argparse, pd, np, st
from utils import ROOT, log_script, log_fn, single_map
from utils.vnt import dtypes, all_vnt_dbs, EFF

# Variant Databases
parser = argparse.ArgumentParser()
parser.add_argument("-H", "--ignore-HGMD", action="store_true", help="Do not process/analyze HGMD Splice variants")
args = parser.parse_args()

vnt_dbs = all_vnt_dbs
if args.ignore_HGMD:
    vnt_dbs.remove("HGMD_splice")

# get NAG strength
strength = pd.read_csv(f"{ROOT}/stats/1nag_motif_freq.txt", sep="\t", index_col=0)["prop_ss"].to_dict()

# load variant events
events = {db : pd.read_csv(f"{ROOT}/variants/found/{db}_nagnag_vnt_events.txt", sep="\t", index_col=0, dtype=dtypes) for db in vnt_dbs}

def diff_strength(x:pd.Series) -> float:

    # alter: alternate vs. reference
    if x["vnt_effect"] == "ALTER":

        r1, r2 = strength[x["ref_seq"][:3].upper()], strength[x["ref_seq"][3:].upper()]
        a1, a2 = strength[x["alt_seq"][:3].upper()], strength[x["alt_seq"][3:].upper()]

        if x["vnt_pos"] < 3:
            return a1 - r1
        return a2 - r2
    
    # create, destroy: affected vs. unaffected
    if x["vnt_effect"] == "CREATE":
        seq_col = "alt_seq"
    elif x["vnt_effect"] == "DESTROY":
        seq_col = "ref_seq"

    s1 = strength[x[seq_col][:3].upper()]
    s2 = strength[x[seq_col][3:].upper()]

    if x["vnt_pos"] < 3:
        return s1 - s2
    return s2 - s1

# NAG strength difference
for db,df in events.items():
    df["dStrength"] = df.apply(lambda x: diff_strength(x), axis=1)
    df.to_csv(f"{ROOT}/variants/found/{db}_nagnag_vnt_events.txt", sep="\t", index=True)
