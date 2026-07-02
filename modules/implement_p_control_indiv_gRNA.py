import numpy as np
import pandas as pd

def implement_p_control_indiv_gRNA(T_gRNA_LFC, binn, p_contr):
    """
    Implement p-controls by LFC targets (instead of p-targets) at gRNA level.
    
    Parameters:
    - T_gRNA_LFC: pd.DataFrame with columns ['gRNA', 'gene', 'LFCM', 'LFCF']
    - binn: array-like, bin edges
    - p_contr: array-like, p/q values corresponding to binn for LFC
    
    Returns:
    - T_q_vertical: pd.DataFrame with columns ['gRNA', 'gene', 'lfc', 'Q']
    """
    #implement p-controls by LFC targets (instead of p-tagets)'

    LFC = T_gRNA_LFC.iloc[:, 2].to_numpy()

    qq = np.zeros(len(LFC))

    for i in range(len(binn) - 1):
        bin_start = binn[i]
        bin_end = binn[i + 1]
        
        # For LFC_1
        idx = (LFC > bin_start) & (LFC <= bin_end)
        qq[idx] = p_contr[i]

    # Assemble new DataFrame
    T_q_vertical = pd.DataFrame({
        'gRNA': T_gRNA_LFC.iloc[:, 0],
        'gene': T_gRNA_LFC.iloc[:, 1],
        'lfc': LFC,
        'Q': qq,
    })

    return T_q_vertical

