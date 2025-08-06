''' Create NAGNAG motif heatmaps, split by splice type '''

from utils import ROOT
from utils.lib import pd, np, plt

def create_motif_heatmap(sfx:str, cm:str, cth:str, tmax:float, nticks:int, cmax:float=None):

    ''' Create 4x4 heatmap of NAGNAG motifs
    
    **Source:** `stats/nagnag_motif_freq.txt`
    
    `figures/heatmaps/{col}_motifs.svg`
    -----------------------------------
    * **x-axis**: identity of N1 (A, C, G, T)
    * **y-axis**: identity of N2 (A, C, G, T)
    * **labels**: number of NAGNAG splice sites of that motif

    Parameters
    ----------
    col : str {"all", "ps", "ds", "as"}
        column of counts to use
    cm : str
        name of Matplotlib colormap to use
    cth : str
        proportion theshold at which labels should be white instead of black
    tmax : int
        maximum tick value to add to colorbar
    nticks : int
        number of ticks to add to colorbar
    cmax : float (default = None)
        maximum proportion to show on colorbar
    '''

    # load motif frequencies
    data = pd.read_csv(f"{ROOT}/stats/nagnag_motif_freq.txt", sep="\t", index_col=0)

    num = data[f"num_{sfx}"].values[:-1].reshape((4,4))
    prop = data[f"prop_{sfx}"].values[:-1].reshape((4,4))

    # figure
    fig = plt.figure(figsize=(12,10), layout="tight")
    ax = fig.gca()

    # heatmap
    im = ax.imshow(prop, cmap=cm, vmin=0, vmax=cmax)

    # annotations
    for i in range(4):
        for j in range(4):
            color = "white" if (prop[i][j] >= cth) else "black"            
            plt.annotate(f"{num[i][j]:,}", xy=(j,i), ha="center", va="center", color=color, size=35)
   
    # colorbar
    cb = fig.colorbar(im, ax=ax, shrink=1, aspect=15, location="right", orientation="vertical", ticks=np.linspace(0, tmax, nticks))
    cb.set_label("Proportion", weight="bold", size=35, labelpad=25)
    cb.ax.tick_params(labelsize=30, pad=10)

    # axis labels
    ax.set_xlabel("N1", weight="bold", size=35, labelpad=15)
    ax.set_ylabel("N2", weight="bold", size=35, labelpad=15)

    # ticks
    ax.tick_params(length=0, pad=10)

    # tick labels
    ax.set_xticks([0, 1, 2, 3], ["A", "C", "G", "T"], size=30)
    ax.set_yticks([0, 1, 2, 3], ["A", "C", "G", "T"], size=30)

    # save figure
    fig.savefig(f"figures/heatmaps/{sfx}_nagnag_motifs.svg", transparent=True)

print("creating motif heatmaps...", end="", flush=True)

create_motif_heatmap(sfx="all", cm="YlOrRd", cth=0.3, tmax=0.3, nticks=4)
create_motif_heatmap(sfx="ps", cm="YlGn", cth=0.35, tmax=0.4, nticks=5)
create_motif_heatmap(sfx="ds", cm="RdPu", cth=0.3, tmax=0.3, nticks=4)
create_motif_heatmap(sfx="as", cm="PuBu", cth=0.35, tmax=0.4, nticks=5)

print("done")
