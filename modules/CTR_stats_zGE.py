import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from .make_histo_crit_stats import make_histo_crit_stats
from .make_LFC_Z_MZ_tables_two import make_LFC_Z_MZ_tables_two
from .q_val_frequentist_critical import q_val_frequentist_critical

def CTR_stats_zGE(
    alf: float,
    st: float,
    en: float,
    step: float,
    T_zGE_12: pd.DataFrame,
    his12t: np.ndarray,
    cond1: str,
    cond2: str,
    condz: str,
    output_dir
):
    """
    Analyse zGE controls vs. targets for two replicates.

    Returns
    -------
    crit_LR12   : (2, ?)  critical LFC limits (rep1, rep2)
    me_sd12     : (2, 2)  mean & SD per replicate
    med_mad12   : (2, 2)  median & MAD per replicate
    binn        : 1‑D bin edges
    p_cont12    : (#bins, 2) control p‑curves (rep1, rep2)
    hiss_cont12 : (#bins, 2) control histograms (rep1, rep2)
    p_targ12    : (#bins, 2) target‑gene p‑curves (rep1, rep2)
    T_zGE_1/2   : DataFrames with LFC, Z, MZ for zGE (rep1 / rep2)
    """

    # ------------------------------------------------------------------
    # 1. Histogram & critical stats for zGE controls
    # ------------------------------------------------------------------
    (
        binn,
        hiss_z12,
        perc_z12,
        crit_1_2_z12,
        bin_p12bz,
        med_mad_z12,
        me_sd_z12,
        mod_z12,
        MZ12,
        Z12,
        MZ_12,
        Z_12,
        n1_n2_n12_z,
    ) = make_histo_crit_stats(
        alf, st, en, step, T_zGE_12, cond1, cond2, output_dir
    )

    me_sd12   = me_sd_z12[:2, :]      # rows 0 & 1  (rep1, rep2)
    med_mad12 = med_mad_z12[:2, :]
    crit_LR12 = crit_1_2_z12[:2, :]

    # ------------------------------------------------------------------
    # 2. zGE tables with LFC / Z / MZ
    # ------------------------------------------------------------------
    T_zGE_1, T_zGE_2 = make_LFC_Z_MZ_tables_two(T_zGE_12, MZ12, Z12)

    # ------------------------------------------------------------------
    # 3. Control p‑curves & histograms (rep1, rep2)
    # ------------------------------------------------------------------
    p_cont12    = bin_p12bz[:, 1:3]    # columns 1 & 2 => rep1 / rep2
    p_cont1, p_cont2 = p_cont12.T
    hiss_cont12 = hiss_z12[:, :2]

    # ------------------------------------------------------------------
    # 4. Target‑gene p‑curves for each replicate
    # ------------------------------------------------------------------
    p_targ1, cL1, cR1, bin_pi1, med_LFCp1, his4p1 = q_val_frequentist_critical(
        alf, binn, his12t[:, 0]
    )
    p_targ2, cL2, cR2, bin_pi2, med_LFCp2, his4p2 = q_val_frequentist_critical(
        alf, binn, his12t[:, 1]
    )
    p_targ12 = np.column_stack([p_targ1, p_targ2])

    # ------------------------------------------------------------------
    # 5. Diagnostic plots
    # ------------------------------------------------------------------
    plt.figure()
    plt.plot(binn, p_cont1, "k--", label="control")
    plt.plot(binn, p_targ1, "k",  label="target")
    plt.title(f"p‑control (--), p‑target (k), indiv gRNA, {cond1}")
    plt.xlabel("LFC bin")
    plt.ylabel("probability to be extreme")
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f"p_distri_targ_cont_{cond1}"), dpi=300, bbox_inches="tight")
    plt.close()
    
    plt.figure()
    plt.plot(binn, p_cont2, "k--", label="control")
    plt.plot(binn, p_targ2, "k",  label="target")
    plt.title(f"p‑control (--), p‑target (k), indiv gRNA, {cond2}")
    plt.xlabel("LFC bin")
    plt.ylabel("probability to be extreme")
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f"p_distri_targ_cont_{cond2}"), dpi=300, bbox_inches="tight")
    plt.close()

    return (
        crit_LR12,
        me_sd12,
        med_mad12,
        binn,
        p_cont12,
        hiss_cont12,
        p_targ12,
        T_zGE_1,
        T_zGE_2,
    )
