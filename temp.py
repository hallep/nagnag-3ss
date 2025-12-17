import pandas as pd
from statsmodels.stats.proportion import proportions_ztest as ztest

def run_ztest(df:pd.DataFrame, col1:str, col2:str):
    
    # z-test
    df["z"], df["p"] = zip(*[ztest([n,c], [sum(df[col1]), sum(df[col2])]) for n,c in zip(df[col1], df[col2])])
    df["p_adj"] = df["p"] * len(df)

    # significance
    df["nStars"] = ((df["p_adj"] < 0.05).astype(int) + (df["p_adj"] < 0.01).astype(int) + 
                        (df["p_adj"] < 0.001).astype(int) + (df["p_adj"] < 0.0001).astype(int))
    df["sig"] = df["nStars"].apply(lambda x: "".join(["*"]*x) if (x > 0) else "ns")

def print_df(df:pd.DataFrame):
    [print("\t".join(map(str, r.values))) for _,r in df.iterrows()]

def nc_v_cds():

    df = pd.read_csv("proteome/nagnag_poswise_freq_count.txt", sep="\t", index_col=0, header=[0,1,2])

    def get_base(base:str) -> pd.DataFrame:
        def get_st(st:str) -> pd.DataFrame:
            dfS = pd.concat((df["nc", st, base], df["cds", st, base]), axis=1, keys=["nc", "cds"])
            dfS["pNC"] = dfS["nc"] / sum(dfS["nc"])
            dfS["pCDS"] = dfS["cds"] / sum(dfS["cds"])
            
            run_ztest(dfS, "nc", "cds")
            return dfS

        stypes = ["PS", "DS", "AS"]
        return pd.concat([get_st(s) for s in stypes], keys=stypes)

    n1 = get_base("n1")
    n2 = get_base("n2")

    print(n1)
    print(n2)

def aatt_by_phase():

    df = pd.read_csv("proteome/transition_type.txt", sep="\t", index_col=[0,1], header=[0,1])

    def get_phase(dfP:pd.DataFrame) -> pd.DataFrame:
        def get_st(st:str):
            return pd.concat((dfP["exp", f"exp_{st}"], dfP["prop", st]), axis=1, keys=["exp", "obs"])

        stypes = ["ps", "ds", "as"]
        return pd.concat([get_st(s) for s in stypes], keys=stypes)

    df0 = get_phase(df.loc[0].loc[["E", "Q", "K", "*"]])
    df1 = get_phase(df.loc[1].loc[["DID", "NID", "CID", "IDR", "NC", "ET"]])
    df2 = get_phase(df.loc[2].loc[["DID", "NID", "CID", "IDR", "NC", "ET"]])

    print_df(df0)
    print_df(df1)
    print_df(df2)

def aa_outcomes():
    df = pd.read_csv("proteome/aa_outcomes.txt", sep="\t", index_col=[0,1,2], header=[0,1])
    dfIns = df.loc[3,"ins",:].loc[["A", "F", "I", "L", "M", "P", "V", "W", "C", "N", "Q", "S", "T", "Y", "D", "E", "H", "K", "R", "G", "*"]]
    dfDel = df.loc[3,"del",:].loc[["A", "F", "I", "L", "M", "P", "V", "W", "C", "N", "Q", "S", "T", "Y", "D", "E", "H", "K", "R", "G", "*"]]

    def get_cols(df:pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        dfS = df[[("exp", "exp_stoch"), ("exp", "exp_splice"), ("prop", "all")]]
        dfS.columns = ["exp_stoch", "exp_splice", "obs"]

        dfT = df[[("prop", "ps"), ("prop", "ds"), ("prop", "as")]]
        dfT.columns = ["ps", "ds", "as"]

        return dfS, dfT

    dfIns1, dfIns2 = get_cols(dfIns)
    dfDel1, dfDel2 = get_cols(dfDel)

    print_df(dfIns1)
    print_df(dfIns2)

    print_df(dfDel1)
    print_df(dfDel2)

