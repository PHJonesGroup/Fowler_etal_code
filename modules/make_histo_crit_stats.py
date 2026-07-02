import pandas as pd
import numpy as np
from .make_histo_LFC import make_histo_LFC
from .compute_p_critLFC import compute_p_critLFC
from .med_mad_MZNP_2 import med_mad_MZNP_2
from .plot_three_panels import plot_three_panels

def make_histo_crit_stats(alf, st, en, step, T_zGE, cond1, cond2, control, output_dir):
    """
    Histograms, critical LFC thresholds, and robust stats for zGE controls,
    for a single contrast (cond1 vs cond2).

    T_zGE : DataFrame with columns [gRNA, gene, lfc]

    Returns
    -------
    bin_      : 1-D bin left-edges
    his       : counts
    perc      : percentages
    crit_LR   : critical LFC thresholds (left/right)
    bin_p     : [bin, p] array
    med_mad   : [median, MAD]
    me_sd     : [mean, SD]
    mod       : mode
    MZ, Z     : MZ / Z scores (1-D, per gRNA)
    n         : count
    """

    # 1. Single LFC vector (3rd column)
    LFC = T_zGE.iloc[:, 2].astype(float).values

    # 2. Histogram
    bin_, his, perc = make_histo_LFC(step, LFC, st, en)

    # 3. Critical LFC thresholds & p-curve (single series)
    bin_p, crit_LR = compute_p_critLFC(alf, bin_, his, cond1, cond2, output_dir)

    # 4. Robust stats, Z / MZ
    med_mad, MZZ, MZ, Z, me_sd, mod = med_mad_MZNP_2(LFC)

    # 5. Diagnostic plot: LFC / Z / MZ distributions
    _, _, perc_z  = make_histo_LFC(step, Z,  st, en)
    _, _, perc_mz = make_histo_LFC(step, MZ, st, en)

    plot_three_panels(
        bin_, perc, perc_z, perc_mz,
        title=f'Distribution of {control} genes: {cond1} vs {cond2}',
        color='b', xlab='LFC',
        save_name=f'distri_{control}_{cond1}_vs_{cond2}', output_dir=output_dir)

    return bin_, his, perc, crit_LR, bin_p, med_mad, me_sd, mod, MZ, Z, len(LFC)