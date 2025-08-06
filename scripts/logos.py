''' Create splice site sequence logos '''

from utils import ROOT
from utils.lib import pd, np, seqlogo
from utils.seq import N, hg38

def get_ss_ppm(sites:pd.DataFrame, site_len:int=0, up_flank:int=0, down_flank:int=0,
               site_col:str="ssite_seq", up_col:str=None, down_col:str=None) -> np.ndarray:

    ''' Get position probability matrix for splice sites
    
    Parameters
    ----------
    sites : pandas.DataFrame
        splice sites/scenarios to consider
    site_len : int (default = 0)
        length of splice site (e.g., 3 for 1-NAG, 6 for NAGNAG)
    up_flank, down_flank : int (default = 0)
        number of bases upstream and downstream (respectively) of splice site to include
    site_col : str (default = "ssite_seq")
        name of column containing splice site sequence
    up_col, down_col : str (default = None)
        name of columns containing upstream and downstream flank sequences (respectively)
    '''

    # get sequences
    up_seq = sites[up_col].str.slice(start=-up_flank) if (up_flank > 0) else pd.Series("", sites.index)
    site_seq = sites[site_col] if (site_len > 0) else pd.Series("", sites.index)
    down_seq = sites[down_col].str.slice(stop=down_flank) if (down_flank > 0) else pd.Series("", sites.index)
    sequence = (up_seq + site_seq + down_seq).str.upper()

    # get poswise base frequencies
    pos = [0] * (up_flank+site_len+down_flank)
    for x in range(len(pos)):
        d = {n : 0 for n in N}
        d.update(sequence.str.get(x).value_counts().to_dict())
        pos[x] = list(d.values())

    # pfm > ppm
    pfm = np.stack(pos)
    return pfm / np.sum(pfm, axis=1, keepdims=True)

# 3' Splice Site Diagrams
def create_3ss_diagrams(iflank:int=30, eflank:int=9):

    ''' Create diagrams for 1-NAG and NAGNAG 3' splice sites
    
    **Source:**
    * `sites/3ss.txt`
    * `sites/nagnag_3ss.txt`

    **Destination:**
    * `figures/logos/1nag_3ss.svg`
    * `figures/logos/nagnag_3ss.svg`
    '''

    # Canonical 1-NAGs
    nag = pd.read_csv(f"{ROOT}/sites/3ss.txt", sep="\t")
    nag = nag[(nag["ssite_type"] == "1C") & (nag["uppercase"] == 1)]

    ppm = get_ss_ppm(nag, site_len=3, up_flank=iflank, down_flank=eflank+3,
                        site_col="ssite_seq", up_col="ssite_iflank_3ss", down_col="ssite_eflank_3ss")
    ppm = np.concatenate((np.array([[1, 0, 0, 0]]), ppm))
    seqlogo.seqlogo(seqlogo.Ppm(ppm), size="xlarge", format="svg", filename=f"{ROOT}/figures/logos/1nag_3ss.svg",
                    show_xaxis=False, show_yaxis=False, stacks_per_line=iflank+eflank+12)

    # NAGNAGs
    nagnag = pd.read_csv(f"{ROOT}/sites/nagnag_3ss.txt", sep="\t")

    ppm = get_ss_ppm(nagnag, site_len=6, up_flank=iflank, down_flank=eflank,
                        site_col="ssite_seq", up_col="ssite_iflank_3ss", down_col="ssite_eflank_3ss")
    ppm = np.concatenate((np.array([[1, 0, 0, 0]]), ppm))
    seqlogo.seqlogo(seqlogo.Ppm(ppm), size="xlarge", format="svg", filename=f"{ROOT}/figures/logos/nagnag_3ss.svg",
                    show_xaxis=False, show_yaxis=False, stacks_per_line=iflank+eflank+12)

# 5' Splice Site Diagrams
def create_5ss_diagrams(eflank:int=5, iflank:int=9, u1_len:int=13):

    ''' Create diagrams for 5' splice sites and the U1 snRNA
    
    **Source:**
    * `sites/3ss.txt`
    * `src/hg38`

    **Destination:**
    * `figures/logos/5ss.svg`
    * `figures/logos/U1snRNA.svg`
    '''
    
    # 5' splice site
    ss = pd.read_csv(f"{ROOT}/sites/3ss.txt", sep="\t")
    ss = ss[ss["uppercase"] == 1]
    cs_inds = [int(i) for inds in ss["csite_inds"] for i in inds.split(",")]
    cs = pd.read_csv(f"{ROOT}/sites/3cs.txt", sep="\t").loc[cs_inds]
    ss5 = pd.DataFrame({
        "eflank" : [e for eflank in cs["scen_eflank_5ss"] for e in eflank.split(",")],
        "iflank" : [i for iflank in cs["scen_iflank_5ss"] for i in iflank.split(",")],
    })

    ppm = get_ss_ppm(ss5, up_flank=eflank, down_flank=iflank, up_col="eflank", down_col="iflank")
    seqlogo.seqlogo(seqlogo.Ppm(ppm), size="xlarge", format="svg", filename=f"{ROOT}/figures/logos/5ss.svg",
                    show_xaxis=False, show_yaxis=False)

    # U1 snRNA
    rnu1 = str(hg38["chr1"][16514121:16514285].complement()).upper()
    u1 = rnu1[-u1_len:]

    ppm = np.zeros((u1_len, 4))
    for x,n in enumerate(u1):
        ppm[x][N.index(n)] = 1

    seqlogo.seqlogo(seqlogo.Ppm(ppm), size="xlarge", format="svg", filename=f"{ROOT}/figures/logos/U1snRNA.svg",
                    ic_scale=False, show_xaxis=False, show_yaxis=False, stack_aspect_ratio=5)

# NAGNAG Logos by Splice Type
def create_nagnag_stype_logo(stype:str, iflank:int=30, eflank:int=6):

    ''' Create motif diagram and site data logo for NAGNAGs of a specific splice type
    
    **Source:** `sites/nagnag_3ss.txt`

    **Destination:**
    * `figures/logos/{stype}_motifs.svg`
    * `figures/logos/{stype}_nagnags.svg`
    '''

    # load sites
    nagnag = pd.read_csv(f"{ROOT}/sites/nagnag_3ss.txt", sep="\t")

    # isolate splice type
    pos = {
        "ps" : "0",
        "ds" : "1",
        "as" : "0,1"
    }[stype]
    sites = nagnag[nagnag["csite_pos"] == pos]

    # motif_logo
    ppm = get_ss_ppm(sites, site_len=6, site_col="ssite_seq")
    seqlogo.seqlogo(seqlogo.Ppm(ppm), size="xlarge", format="svg", filename=f"{ROOT}/figures/logos/{stype}_motifs.svg",
                    ic_scale=False, show_xaxis=False, show_yaxis=False)

    # motif_logo
    ppm = get_ss_ppm(sites, site_len=6, up_flank=iflank, down_flank=eflank,
                     site_col="ssite_seq", up_col="ssite_iflank_3ss", down_col="ssite_eflank_3ss")
    seqlogo.seqlogo(seqlogo.Ppm(ppm), size="xlarge", format="svg", filename=f"{ROOT}/figures/logos/{stype}_nagnags.svg",
                    first_index=-(iflank+6))

print("creating sequence logos...", end="", flush=True)

create_3ss_diagrams()
create_5ss_diagrams()

create_nagnag_stype_logo("ps")
create_nagnag_stype_logo("ds")
create_nagnag_stype_logo("as")

print("done")
