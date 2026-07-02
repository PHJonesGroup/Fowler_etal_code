import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from .make_histo_crit_stats import make_histo_crit_stats
from .make_LFC_Z_MZ_tables_two import make_LFC_Z_MZ_tables_two
from .q_val_frequentist_critical import q_val_frequentist_critical

def CTR_stats_zGE(alf, st, en, step, T_zGE, hist, cond1, cond2, condz, control, output_dir):
    """
    Analyse zGE controls vs. targets for a single contrast (cond1 vs cond2).

    Returns
    -------
    crit_LR   : critical LFC limits (left/right)
    me_sd     : (mean, SD) of the zGE LFC
    med_mad   : (median, MAD) of the zGE LFC
    binn      : bin edges
    p_cont    : control p-curve
    hiss_cont : control histogram
    p_targ    : target-gene p-curve
    T_zGE_tab : zGE table with LFC, Z, MZ
    """
    print(T_zGE)
    # 1. Histogram & critical stats for zGE controls
    (binn, hiss_z, perc, crit_LR, bin_pz, med_mad_z, me_sd_z, mod, MZ, Z, n_z
    ) = make_histo_crit_stats(alf, st, en, step, T_zGE, cond1, cond2, control, output_dir)

    me_sd   = me_sd_z
    med_mad = med_mad_z
    crit_LR = crit_LR
    # 2. zGE table with LFC / Z / MZ
    T_zGE_tab = make_LFC_Z_MZ_tables_two(T_zGE, MZ, Z)

    # 3. Control p-curve & histogram
    p_cont    = bin_pz[:, 1]     # single control p-curve column
    hiss_cont = hiss_z

    # 4. Target-gene p-curve
    p_targ, cL, cR, bin_pi, med_LFCp, his4p = q_val_frequentist_critical(alf, binn, hist)

    # 5. Plot: target vs. control p-curves
    plt.figure()
    plt.plot(binn, p_cont, color="0.5", linestyle="--", label=f"control {control}")
    plt.plot(binn, p_targ, color="C3", linestyle="-",  label="target")
    plt.title(f"Target vs. {control} gRNAs ({cond1} vs {cond2})")
    plt.xlabel("LFC bin")
    plt.ylabel("Probability")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f"p_distri_targ_cont_{cond1}_vs_{cond2}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    return crit_LR, me_sd, med_mad, binn, p_cont, hiss_cont, p_targ, T_zGE_tab