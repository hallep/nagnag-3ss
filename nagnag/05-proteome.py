''' Calculate amino acid transitions for NAGNAG splice scenarios

`sites/nagnag_scens.txt`
------------------------
* **vc_ps** (*str*): variable codons when proximally spliced
* **vc_ds** (*str*): variable codons when distally spliced
* **aa_ps** (*str*): variable amino acid(s) when proximally spliced
* **aa_ds** (*str*): variable amino acid when distally spliced
* **aattype** (*str*): amino acid transition type
* **ins_aa** (*str*): inserted amino acids (from distal to proximal)
* **del_aa** (*str*): deleted amino acids (from distal to proximal)
'''

from lib import pd
from utils import ROOT, log_script, log_fn
from utils.seq import get_tsn, get_categorical_outcomes

# ===== RUN ===== #
log_script("05-proteome.py")
log_fn("Calculating proteomic effects")

# load splice scenarios
scens = pd.read_csv(f"{ROOT}/sites/nagnag_scens.txt", sep="\t", index_col=0)

# get transitions
tsn = [get_tsn(u,m,d,p) for u,m,d,p in zip(scens["up_seq"], scens["ssite_seq"],
                                                scens["down_seq"], scens["phase"])]
vc_ps, vc_ds, aa_ps, aa_ds, aattype = zip(*tsn)

ins_aa, del_aa = zip(*[get_categorical_outcomes(p, ps, ds) for p,ps,ds in zip(scens["phase"], aa_ps, aa_ds)])

# add columns to DataFrame
scens["vc_ps"] = vc_ps
scens["vc_ds"] = vc_ds
scens["aa_ps"] = aa_ps
scens["aa_ds"] = aa_ds
scens["aattype"] = aattype
scens["ins_aa"] = ins_aa
scens["del_aa"] = del_aa

# save DataFrame
scens.to_csv(f"{ROOT}/sites/nagnag_scens.txt", sep="\t", index=True)
