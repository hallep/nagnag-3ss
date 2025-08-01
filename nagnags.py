''' Run all NAGNAG analyses '''

import os

# Get data
def data():

    ''' Get transcript annotation and chromosome sequence data '''

    import scripts.data as data

    # make "source" folder
    os.makedirs("source", exist_ok=True)

    # download NCBI RefSeq data
    data.download_ncbi_refseq()

    # download hg38 sequences
    data.download_hg38()

    # create pickle of sequence dictionary
    data.create_hg38_seq_dict()

# Parse transcript annotations
def tx_annotations():

    ''' Parse NCBI RefSeq transcript annotation data '''

    import scripts.refseq as refseq

    # make "data" folder
    os.makedirs("data", exist_ok=True)

    # extract cleavage sites
    refseq.get_cleavage_sites()

    # add cleavage site flanking sequences
    refseq.add_csite_sequences()

    # group cleavage sites into splice sites
    refseq.group_csites()

    # remove splice sites in repeat masker
    refseq.identify_uppercase()

    # get NAGNAGs
    refseq.extract_nagnags()

    # get NAGNAG splice scenarios
    refseq.split_nagnag_scenarios()

    # add splice scenario flanking sequences
    refseq.add_scen_flanks()

# Get site frequencies
def statistics():

    ''' Get splice site, motif, splice type, and phase frequencies '''

    import scripts.stats as stats

    # make "stats" folder
    os.makedirs("stats", exist_ok=True)

    # splice site breakdown (Table 1)
    stats.count_ss_types()

    # canonical NAG motifs (Table 4)
    stats.count_1nag_motifs()

    # NAGNAG splice types (Table 2)
    stats.count_nagnag_splice_types()

    # NAGNAG splice scenario phases (Table 3)
    stats.count_nagnag_phases()

# Analyze NAGNAG motifs by splice site
def motifs():

    ''' Analyze NAGNAG motifs '''

    import scripts.motifs as motifs

    # make "stats" folder
    os.makedirs("stats", exist_ok=True)

    # make "figures" folder
    os.makedirs("figures", exist_ok=True)

    # count NAGNAG motifs (Tables S1-3)
    motifs.count_motifs()

    # create motif heatmaps (Figure 4A, C, E)
    motifs.create_motif_heatmap(col="all", cm="YlOrRd", cth=0.3, tmax=0.3, nticks=4)
    motifs.create_motif_heatmap(col="ps", cm="YlGn", cth=0.35, tmax=0.4, nticks=5)
    motifs.create_motif_heatmap(col="ds", cm="RdPu", cth=0.3, tmax=0.3, nticks=4)
    motifs.create_motif_heatmap(col="as", cm="PuBu", cth=0.35, tmax=0.4, nticks=5)

# Create splice site logos
def logos():

    ''' Create splice site sequence logos '''

    import scripts.logos as logos

    # make "figures" folder
    os.makedirs("figures", exist_ok=True)

    # create graphical logo for 1-NAG 3' splice sites (Figure 1A)
    logos.create_1nag_logo()

    # create graphical logo for NAGNAG 3' splice sites (Figure 1B)
    logos.create_nagnag_logo()

    # create logos for each splice type (Figure 4B, D, F)
    logos.create_nagnag_splice_type_logos()

# Analyze proteome effects
def proteome():

    ''' Analyze the effects of NAGNAG alternative splicing on the proteome '''

    import scripts.aa as aa

    # make "protein" folder
    os.makedirs("protein", exist_ok=True)

    # determine NAGNAG amino acid transitions
    aa.add_nagnag_aa_transitions()

    # define all possible variable codons
    aa.get_possible_vc()

    # count NAGNAG variable codons
    aa.count_nagnag_vc(0)
    aa.count_nagnag_vc(1)
    aa.count_nagnag_vc(2)

    # count NAGNAG amino acid transitions (Tables S4-6)
    aa.count_aa_transitions(0)
    aa.count_aa_transitions(1)
    aa.count_aa_transitions(2)

    # group by amino acid transition types
    aa.add_aatype(1)
    aa.add_aatype(2)

    # count amino acid transition types (Tables 5 and 6)
    aa.count_aatype(1)
    aa.count_aatype(2)

data()
tx_annotations()
statistics()
motifs()
logos()
proteome()
