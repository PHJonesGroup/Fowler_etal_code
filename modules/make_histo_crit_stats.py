import pandas as pd
import numpy as np
from .make_histo_vec_rep12 import make_histo_vec_rep12
from .compute_p_critLFC import compute_p_critLFC
from .med_mad_MZNP_2 import med_mad_MZNP_2
from .plot_three_panels import plot_three_panels

def make_histo_crit_stats(
    alf: float,
    st: float,
    en: float,
    step: float,
    T_zGE_12: pd.DataFrame,
    cond1: str,
    cond2: str,
    output_dir
):
    """
    Compute histograms, critical LFC thresholds, and robust stats for zGE
    controls in two replicates plus the pooled set.

    Parameters
    ----------
    alf   : significance level (e.g. 0.06)
    st/en : histogram lower / upper bounds for LFC
    step  : histogram bin width
    T_zGE_12 : DataFrame with columns [gRNA, gene, lfc1, lfc2]
    cond1, cond2 : str
        Condition label for plot title.
    Returns 
    -------
    bin           : 1‑D ndarray of bin left‑edges
    his           : (#bins, 3) counts for rep1, rep2, pooled
    perc          : (#bins, 3) percentages for rep1, rep2, pooled
    crit_1_2_12   : (3, ?) matrix of critical LFC thresholds
    bin_p12b      : (#bins, 4)  [bin  p_rep1  p_rep2  p_pooled]
    med_mad       : (3, 2)   [median, MAD] for rep1, rep2, pooled
    me_sd         : (3, 2)   [mean,   SD]  for rep1, rep2, pooled
    mod           : length‑3 array of modes
    MZ, Z         : (#rows, 2) MZ / Z for rep1 & rep2
    MZ_12, Z_12   : 1‑D arrays for pooled replicate
    n1_n2_n12     : (#bins, 3) cumulative counts (from make_histo_vec_rep12)
    """

    # ------------------------------------------------------------------
    # --- 1.  Pull LFC vectors for two replicates -----------------------
    # ------------------------------------------------------------------
    LFC_1 = T_zGE_12.iloc[:, 2].astype(float).values
    LFC_2 = T_zGE_12.iloc[:, 3].astype(float).values
    LFC_12 = np.concatenate([LFC_1, LFC_2])

    # ------------------------------------------------------------------
    # --- 2.  Histograms for rep1 / rep2 / pooled ----------------------
    # ------------------------------------------------------------------
    (
        bin_edges,
        perc_1,
        perc_2,
        perc_12,
        his_1,
        his_2,
        his_12,
        n1_n2_n12,
    ) = make_histo_vec_rep12(step, LFC_1, LFC_2, st, en)

    perc = np.column_stack([perc_1, perc_2, perc_12])
    his  = np.column_stack([his_1,  his_2,  his_12])
    bin_ = bin_edges 

    # ------------------------------------------------------------------
    # --- 3.  Critical LFC thresholds & p‑curves ------------------------
    # ------------------------------------------------------------------
    bin_p12b, crit_1_2_12 = compute_p_critLFC(
        alf, bin_, his_1, his_2, his_12, cond1, cond2, output_dir)

    # ------------------------------------------------------------------
    # --- 4.  Med‑MAD, Mean‑SD, Mode, Z‑scores, MZ‑scores --------------
    # ------------------------------------------------------------------
    med_mad_1, MZZ1, MZ_1, Z_1, me_sd_1, mod_1 = med_mad_MZNP_2(LFC_1)
    med_mad_2, MZZ2, MZ_2, Z_2, me_sd_2, mod_2 = med_mad_MZNP_2(LFC_2)
    med_mad_12, MZZ12, MZ_12, Z_12, me_sd_12, mod_12 = med_mad_MZNP_2(LFC_12)

    med_mad = np.vstack([med_mad_1, med_mad_2, med_mad_12])
    me_sd   = np.vstack([me_sd_1,   me_sd_2,   me_sd_12])
    mod     = np.array([mod_1, mod_2, mod_12])

    # Stack MZ/Z for two reps (pooled goes out separately)
    MZ = np.column_stack([MZ_1, MZ_2])
    Z  = np.column_stack([Z_1,  Z_2])

    # ------------------------------------------------------------------
    # --- 5.  Diagnostic plots ---------------------------
    # ------------------------------------------------------------------
    # --- Compute Z and MZ histograms ---
    _, perc_z1, perc_z2, perc_z12, _, _, _, _ = make_histo_vec_rep12(step, Z_1, Z_2, st, en)
    _, perc_mz1, perc_mz2, perc_mz12, _, _, _, _ = make_histo_vec_rep12(step, MZ_1, MZ_2, st, en)

    # LFC / Z / MZ
    # replicate‑1
    plot_three_panels(bin_, perc_1, perc_z1, perc_mz1,
                   title=f'Distribution of zero expressed genes in condition: {cond1}', color='b', xlab='LFC', save_name=f'distri_zGE_{cond1}', output_dir=output_dir)

    # replicate‑2
    plot_three_panels(bin_, perc_2, perc_z2, perc_mz2,
                   title=f'Distribution of zero expressed genes in condition: {cond2}', color='r', xlab='LFC', save_name=f'distri_zGE_{cond2}', output_dir=output_dir)

    # ------------------------------------------------------------------
    return (
        bin_,
        his,
        perc,
        crit_1_2_12,
        bin_p12b,
        med_mad,
        me_sd,
        mod,
        MZ,
        Z,
        MZ_12,
        Z_12,
        n1_n2_n12,
    )