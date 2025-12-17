import pandas as pd
from statsmodels.stats.proportion import proportions_ztest as ztest

df = pd.read_csv("proteome/nagnag_poswise_freq_count.txt", sep="\t", index_col=0, header=[0,1,2])

def get_base(base:str) -> pd.DataFrame:
    def get_st(st:str) -> pd.DataFrame:
        dfS = pd.concat((df["nc", st, base], df["cds", st, base]), axis=1, keys=["nc", "cds"])
        dfS["pNC"] = dfS["nc"] / sum(dfS["nc"])
        dfS["pCDS"] = dfS["cds"] / sum(dfS["cds"])

        # z-test
        dfS["z"], dfS["p"] = zip(*[ztest([n,c], [sum(dfS["nc"]), sum(dfS["cds"])]) for n,c in zip(dfS["nc"], dfS["cds"])])
        dfS["p_adj"] = dfS["p"] * 4

        dfS["nStars"] = ((dfS["p_adj"] < 0.05).astype(int) + (dfS["p_adj"] < 0.01).astype(int) + 
                         (dfS["p_adj"] < 0.001).astype(int) + (dfS["p_adj"] < 0.0001).astype(int))
        dfS["sig"] = dfS["nStars"].apply(lambda x: "".join(["*"]*x) if (x > 0) else "ns")

        return dfS

    stypes = ["PS", "DS", "AS"]
    return pd.concat([get_st(s) for s in stypes], keys=stypes)

n1 = get_base("n1")
n2 = get_base("n2")

print(n1)
