import pandas as pd

def make_tables_Z_two(target_12, Z_t1, Z_t2, cond1, cond2):
    """
    Takes a vertical gRNA gene table with LFC1 & LFC2, Q1 & Q2,
    adds Z1, Z2 columns, and makes two separate tables with Z for rep1 and rep2.
    
    Parameters:
    - target_12: pd.DataFrame with columns [gRNA, gene, lfc_1, lfc_2, Q1, Q2]
    - Z_t1, Z_t2: numpy arrays or lists of Z scores for replicate 1 and 2
    - cond1, cond2: string labels for conditions for rep1 and rep2
    
    Returns:
    - T_t1_wt, T_t2_wt: pd.DataFrames with columns [gRNA, gene, LFC, Z_zGE_LFC, Q, condition]
    """

    # Extract columns
    gRNA = target_12.iloc[:, 0]
    gene = target_12.iloc[:, 1]
    LFC_t1 = target_12.iloc[:, 2]
    LFC_t2 = target_12.iloc[:, 3]
    Q1 = target_12.iloc[:, 4]
    Q2 = target_12.iloc[:, 5]

    # Replicate 1 table
    condition1 = [cond1] * len(LFC_t1)
    T_t1_wt = pd.DataFrame({
        'gRNA': gRNA,
        'gene': gene,
        'LFC': LFC_t1,
        'Z_zGE_LFC': Z_t1,
        'Q': Q1,
        'condition': condition1
    })

    # Replicate 2 table
    condition2 = [cond2] * len(LFC_t2)
    T_t2_wt = pd.DataFrame({
        'gRNA': gRNA,
        'gene': gene,
        'LFC': LFC_t2,
        'Z_zGE_LFC': Z_t2,
        'Q': Q2,
        'condition': condition2
    })

    return T_t1_wt, T_t2_wt
