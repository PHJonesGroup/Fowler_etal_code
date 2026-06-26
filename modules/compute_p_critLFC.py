import numpy as np
import os
import matplotlib.pyplot as plt
from .q_val_frequentist_critical import q_val_frequentist_critical

def compute_p_critLFC(alf, binn, hiss_1, hiss_2, hiss_12, cond1, cond2, output_dir):
    """
    Compute p-values and critical values for LFC distributions of replicates and pooled.

    Parameters
    ----------
    alf : float
        Significance level.
    binn : array_like
        Bin edges or centers for LFC histogram.
    hiss_1, hiss_2, hiss_12 : array_like
        Histogram counts for replicate 1, replicate 2, and pooled data.
    cond1, cond2 : str
        Condition label for plot title.

    Returns
    -------
    bin_pzit : ndarray
        Array with columns: [bin, p_rep1, p_rep2, p_pooled]
    crit_12_targ : ndarray
        Critical values for rep1, rep2, and pooled (shape 3x2)
    """

    # Replicate 1
    p1, clz, crz, bin_pz, med_LFCpz, hiss4pz = q_val_frequentist_critical(alf, binn, hiss_1)
    critz = [clz, crz]

    # Replicate 2 (intergenic)
    p2, cli, cri, bin_pi, med_LFCpi, hiss4pi = q_val_frequentist_critical(alf, binn, hiss_2)
    criti = [cli, cri]

    # Both replicates pooled
    p12, clT, crT, bin_pT, med_LFCpT, hiss4pT = q_val_frequentist_critical(alf, binn, hiss_12)
    critT = [clT, crT]

    crit_12_targ = np.array([critz, criti, critT])

    # Assemble bin_pzit matrix with bin and p-values for rep1, rep2, pooled
    bin_pzit = np.column_stack((binn, p1, p2, p12))

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.plot(binn, p1, 'b', label={cond1})
    plt.plot(binn, p12, 'k:', label='Both')
    plt.plot(binn, p2, 'r', label={cond2})
    plt.grid(True)
    plt.xlabel('LFC')
    plt.title(f'zGE distribution', fontsize=14)
    plt.legend()
    plt.savefig(os.path.join(output_dir, f"p_controls_zGE_{cond1}{cond2}"), dpi=300, bbox_inches="tight")
    plt.close()

    return bin_pzit, crit_12_targ

