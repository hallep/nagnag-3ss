''' Analyze NAGNAG motifs by splice type '''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scripts.sequence import get_NAGNAGs

# Count NAGNAG motifs
def count_motifs():

    ''' Count the number of each NAGNAG motif (16 in total)
    
    src: data/nagnag_3ss.txt

    dst: stats/motif_count.txt

        "motif" : str
            NAGNAG motif
        "all" : int
            number of NAGNAG splice sites of each motif
        "ps" : int
            number of proximally-spliced NAGNAG splice sites of each motif
        "ds" : int
            number of distally-spliced NAGNAG splice sites of each motif
        "as" : int
            number of alternatively-spliced NAGNAG splice sites of each motif
    '''

    # load splice sites
    ssites = pd.read_csv("data/nagnag_3ss.txt", sep="\t", index_col=0)

    motif = get_NAGNAGs()
    count_all = [len(ssites[ssites["ssite_seq"] == m]) for m in motif]
    count_ps = [len(ssites[(ssites["csite_pos"] == "0") & (ssites["ssite_seq"] == m)]) for m in motif]
    count_ds = [len(ssites[(ssites["csite_pos"] == "1") & (ssites["ssite_seq"] == m)]) for m in motif]
    count_as = [len(ssites[(ssites["csite_pos"] == "0,1") & (ssites["ssite_seq"] == m)]) for m in motif]

    # create DataFrame
    counts = pd.DataFrame({
        "motif" : motif,
        "all" : count_all,
        "ps" : count_ps,
        "ds" : count_ds,
        "as" : count_as
    })

    # save DataFrame
    counts.to_csv("stats/motif_count.txt", sep="\t", index=False)

# Show NAGNAG motifs as heatmap
def create_motif_heatmap(col:str, cm:str, cth:float, nticks:int, tmax:int, tmin:int=0, pmax:float=None, raw=True):
    
    ''' Create 4x4 heatmap of NAGNAG motifs
    
    x-axis: identity of N1 (A, C, G, T)
    y-axis: identity of N2 (A, C, G, T)

    src: stats/motif_count.txt
    
    dst: figures/{col}_motifs_{raw}.svg
    
    Parameters
    ----------
    col : str {"all", "ps", "ds", "as"}
        column of counts to use
    cm : str
        name of Matplotlib colormap to use
    cth : str
        proportion theshold at which labels should be white instead of black
    nticks : int
        number of ticks to add to colorbar
    tmax : int
        maximum tick value to add to colorbar
    tmin : int (default = 0)
        minimum tick value to add to colorbar
    pmax : float (default = None)
        maximum proportion to show on colorbar
    raw : bool (default = True)
        format of value annotations
         * True: show raw counts
         * False: show proportions
    '''

    # annotation (ats), label (lts), and tick (tts) text sizes
    ats = 35
    lts = 35
    tts = 30

    # heatmap axis (hlp), colorbar axis (clp), and tick (tlp) label padding
    hlp = 15
    clp = 25
    tlp = 10

    # load data
    mfreq = pd.read_csv("stats/motif_count.txt", sep="\t")

    # get values
    count = np.array([mfreq[col][0:4], mfreq[col][4:8], mfreq[col][8:12], mfreq[col][12:16]]).transpose()
    prop = count / count.sum()
    
    # figure
    fig = plt.figure(figsize=(12,10), layout="tight")
    ax = fig.gca()

    # heatmap
    im = ax.imshow(prop, cmap=cm, vmin=0, vmax=pmax)

    # annotations
    for i in range(4):
        for j in range(4):

            # define label color
            if prop[i][j] >= cth:
                color = "white"
            else:
                color = "black"
            
            if raw:
                plt.annotate(f"{count[i][j]:,}", xy=(j,i), ha="center", va="center", color=color, size=ats)
            else:
                plt.annotate(f"{prop[i][j]:,}", xy=(j,i), ha="center", va="center", color=color, size=ats)
   
    # colorbar
    cb = fig.colorbar(im, ax=ax, shrink=1, aspect=15, location="right", orientation="vertical", ticks=np.linspace(tmin, tmax, nticks))
    cb.set_label("Proportion", weight="bold", size=lts, labelpad=clp)
    cb.ax.tick_params(labelsize=tts, pad=tlp)

    # axis labels
    ax.set_xlabel("N1", weight="bold", size=lts, labelpad=hlp)
    ax.set_ylabel("N2", weight="bold", size=lts, labelpad=hlp)

    # ticks
    ax.tick_params(length=0, pad=tlp)

    # tick labels
    ax.set_xticks([0, 1, 2, 3], ["A", "C", "G", "T"], size=tts)
    ax.set_yticks([0, 1, 2, 3], ["A", "C", "G", "T"], size=tts)

    plt.show()

    # save figure
    fig.savefig(f"figures/{col}_motifs_{"raw" if raw else "prop"}.svg")
