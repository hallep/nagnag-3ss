''' Create splice site sequence logos '''

import pandas as pd
import numpy as np
import seqlogo

# Create position probability matrix
def create_ppm(pfm:np.ndarray) -> np.ndarray:

    ''' Create a vertical position matrix

    Parameters
    ----------
    pfm : numpy.NDArray
        position frequency matrix
         * counts for each base at every position
         * size: {length} x 4
    
    Returns
    -------
    ppm : numpy.NDArray
        position probability matrix
         * proportions for each base at every position
         * size: {length} x 4
    '''

    return pfm / np.sum(pfm, axis=1, keepdims=True)

# Create position frequency matrix for 3' splice sites
def create_3ss_pfm(sites:pd.DataFrame, slen:int, iflank:int, eflank:int, seq_col="ssite_seq",
                   iseq_col="ssite_iflank_3ss", eseq_col="ssite_eflank_3ss") -> tuple[np.ndarray, np.ndarray]:
    
    ''' Create PFM and PPM for 3' splice sites

    Parameters
    ----------
    sites : pandas.DataFrame
        DataFrame of splice site information
    slen : int
        length of splice site motif
    iflank : int
        number of upstream intronic bases to include
    eflank : int
        number of downstream exonic bases to include
    ss_file : str (default = "/mnt/data_3/hallep/thesis/nagnag_3ss.txt")
        filepath to tab-delimited splice site file
    seq_col : str (default = "ssite_seq")
        name of column containing splice site sequence
    iseq_col : str (default = "ssite_iflank_3ss")
        name of column containing intronic sequence upstream of splice site
    eseq_col : str (default = "ssite_eflank_3ss")
        name of column containing exonic sequence downstream of splice site
    
    Returns
    -------
    pfm : numpy.NDArray
        vertical position frequency matrix
    ppm : numpy.NDArray
        vertical position probability matrix
    '''

    num_bases = iflank + slen + eflank
    pfm = np.array([[0] * 4 for _ in range(num_bases)])

    base_dict = {
        "A":0,
        "C":1,
        "G":2,
        "T":3
    }

    # for each splice site
    for _,ss in sites.iterrows():

        # upstream intronic sequence
        if iflank == 0:
            iseq = ""
        else:
            iseq =  ss[iseq_col][-iflank:]

        # downstream exonic sequence
        if eflank == 0:
            eseq = ""
        else:
            eseq = ss[eseq_col][:eflank]
        
        # sequence of interest
        seq = iseq + ss[seq_col] + eseq

        # for each base in sequence
        for pos,b in enumerate(seq):
            base = base_dict[b.upper()]

            # add counts to PFMs
            pfm[pos][base] += 1

    # get position probability matrix
    ppm = create_ppm(pfm)

    return pfm, ppm

# Create canonical 1-NAG  3' splice site logo
def create_1nag_logo(iflank:int=30, eflank:int=5):

    ''' Create a graphical sequence logo for 1-NAG 3' splice sites, including the upstream intronic branch point
    
    src: data/3ss.txt
    
    dst: figures/1nag_{iflank}u{eflank}d_3ss_logo.svg
    
    Parameters
    ----------
    iflank : int (default = 30)
        number of intronic bases upstream of the NAG to show
    eflank : int (default = 5)
        number of exonic bases downstream of the NAG to show
    '''

    # load data
    ssites = pd.read_csv("data/3ss.txt", sep="\t", index_col=0)
    ssites = ssites[(ssites["ssite_type"] == "1C") & (ssites["uppercase"] == 1)]

    # get position probability matrix
    _, raw_ppm = create_3ss_pfm(ssites, 3, iflank, eflank)

    # add upstream intronic branchpoint
    ppm_bp = np.concatenate((np.array([[1, 0, 0, 0]]), raw_ppm))

    # create Ppm object
    ppm = seqlogo.Ppm(ppm_bp)

    # create logo
    seqlogo.seqlogo(ppm, size="xlarge", format="svg",
                    filename=f"figures/1nag_{iflank}u{eflank}d_3ss_logo.svg",
                    show_xaxis=False, show_yaxis=False)

# Create NAGNAG 3' splice site logo
def create_nagnag_logo(iflank:int=30, eflank:int=5):

    ''' Create a graphical sequence logo for NAGNAG 3' splice sites, including the upstream intronic branch point
    
    src: data/nagnag_3ss.txt
    
    dst: figures/nagnag_{iflank}u{eflank}d_3ss_logo.svg
    
    Parameters
    ----------
    iflank : int (default = 30)
        number of intronic bases upstream of the NAG to show
    eflank : int (default = 5)
        number of exonic bases downstream of the NAG to show
    '''

    # load data
    ssites = pd.read_csv("data/nagnag_3ss.txt", sep="\t", index_col=0)

    # get position probability matrix
    _, raw_ppm = create_3ss_pfm(ssites, 6, iflank, eflank)

    # add upstream intronic branchpoint
    ppm_bp = np.concatenate((np.array([[1, 0, 0, 0]]), raw_ppm))

    # create Ppm object
    ppm = seqlogo.Ppm(ppm_bp)

    # create logo
    seqlogo.seqlogo(ppm, size="xlarge", format="svg",
                    filename=f"figures/nagnag_{iflank}u{eflank}d_3ss_logo.svg",
                    show_xaxis=False, show_yaxis=False, stacks_per_line=45)

# Create NAGNAG sequence logos by splice type
def create_nagnag_splice_type_logos(iflank:int=24, eflank:int=5):

    ''' Create a sequence logos for NAGNAGs by splice type
    
    src: data/nagnag_3ss.txt
    
    dst: figures/nagnag_{iflank}u{eflank}d_3ss_logo.svg
    
    Parameters
    ----------
    iflank : int (default = 30)
        number of intronic bases upstream of the NAG to show
    eflank : int (default = 5)
        number of exonic bases downstream of the NAG to show
    '''

    # load and split data
    ssites = pd.read_csv("/mnt/data_3/hallep/thesis/nagnag_3ss.txt", sep="\t", index_col=0)
    ps_ssites = ssites[ssites["csite_pos"] == "0"]
    ds_ssites = ssites[ssites["csite_pos"] == "1"]
    as_ssites = ssites[ssites["csite_pos"] == "0,1"]

    # create position probability matrices
    ps_ppm = seqlogo.Ppm(create_3ss_pfm(ps_ssites, 6, iflank, eflank)[1])
    ds_ppm = seqlogo.Ppm(create_3ss_pfm(ds_ssites, 6, iflank, eflank)[1])
    as_ppm = seqlogo.Ppm(create_3ss_pfm(as_ssites, 6, iflank, eflank)[1])

    # proximally-spliced NAGNAGs
    seqlogo.seqlogo(ps_ppm, size="xlarge", format="svg",
                    filename=f"figures/nagnag_{iflank}u{eflank}d_ps_logo.svg",
                    first_index=-(iflank+6))

    # proximally-spliced NAGNAGs
    seqlogo.seqlogo(ds_ppm, size="xlarge", format="svg",
                    filename=f"figures/nagnag_{iflank}u{eflank}d_ds_logo.svg",
                    first_index=-(iflank+6))

    # alternatively-spliced NAGNAGs
    seqlogo.seqlogo(as_ppm, size="xlarge", format="svg",
                    filename=f"figures/nagnag_{iflank}u{eflank}d_as_logo.svg",
                    first_index=-(iflank+6))
