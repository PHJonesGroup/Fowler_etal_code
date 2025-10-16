import numpy as np
import pandas as pd

def implement_p_control_indiv_gRNA(T_gRNA_LFC, binn, p_contr1, p_contr2):
    """
    Implement p-controls by LFC targets (instead of p-targets) at gRNA level.
    
    Parameters:
    - T_gRNA_LFC: pd.DataFrame with columns ['gRNA', 'gene', 'LFCM', 'LFCF']
    - binn: array-like, bin edges
    - p_contr1: array-like, p/q values corresponding to binn for LFC_1
    - p_contr2: array-like, p/q values corresponding to binn for LFC_2
    
    Returns:
    - T_12_q_vertical: pd.DataFrame with columns ['gRNA', 'gene', 'lfc_1', 'lfc_2', 'Q1', 'Q2']
    """
    #implement p-controls by LFC targets (instead of p-tagets)'

    LFC_1 = T_gRNA_LFC.iloc[:, 2].to_numpy()  # LFCM
    LFC_2 = T_gRNA_LFC.iloc[:, 3].to_numpy()  # LFCF

    qq1 = np.zeros(len(LFC_1))
    qq2 = np.zeros(len(LFC_2))

    for i in range(len(binn) - 1):
        bin_start = binn[i]
        bin_end = binn[i + 1]
        
        # For LFC_1
        idx_1 = (LFC_1 > bin_start) & (LFC_1 <= bin_end)
        qq1[idx_1] = p_contr1[i]

        # For LFC_2
        idx_2 = (LFC_2 > bin_start) & (LFC_2 <= bin_end)
        qq2[idx_2] = p_contr2[i]

    # Assemble new DataFrame
    T_12_q_vertical = pd.DataFrame({
        'gRNA': T_gRNA_LFC.iloc[:, 0],
        'gene': T_gRNA_LFC.iloc[:, 1],
        'lfc_1': LFC_1,
        'lfc_2': LFC_2,
        'Q1': qq1,
        'Q2': qq2,
    })

    return T_12_q_vertical

