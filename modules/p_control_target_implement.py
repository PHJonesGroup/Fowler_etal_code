import pandas as pd
from .implement_p_control_indiv_gRNA import implement_p_control_indiv_gRNA

def p_control_target_implement(T_target, T_zGE, binn, p_cont12):
    """
    Given p per bin, infer p/q values for [T_target; T_zGE] from LFC values.

    Parameters:
    - T_target: pd.DataFrame with gRNA target data
    - T_zGE: pd.DataFrame with zGE gene data
    - binn: array-like, bin edges
    - p_cont12: 2D array-like, calibration p/q values per bin (columns for rep1 and rep2)

    Returns:
    - T_ij_q_verti: pd.DataFrame with gRNAs, genes, lfc_1, lfc_2, Q1, Q2
    """

    #FDR  FRR (false rejection)-correct p-values for both tails of LFC NNT distribution
    #implement p-controls into LFC (instead of p-targets)

    # Concatenate tables vertically 
    T_target_zGE = pd.concat([T_target, T_zGE], ignore_index=True)

    # Extract p/q control columns
    p_contr1 = p_cont12[:, 0]
    p_contr2 = p_cont12[:, 1]

    # Call the function that assigns p/q values based on LFC bins
    T_ij_q_verti = implement_p_control_indiv_gRNA(T_target_zGE, binn, p_contr1, p_contr2)
    
    return T_ij_q_verti