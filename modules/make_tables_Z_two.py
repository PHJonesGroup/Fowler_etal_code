import pandas as pd

def make_tables_Z_two(target, Z_t):
    """
    Takes a vertical gRNA gene table with LFC, Q1
    adds Z columns, and makes two separate tables with Z .
    
    Parameters:
    - target: pd.DataFrame with columns [gRNA, gene, lfc, Q1]
    - Z_t: numpy arrays or lists of Z scores 
    
    Returns:
    - T_wt: pd.DataFrames with columns [gRNA, gene, LFC, Z_zGE_LFC, Q]
    """

    # Extract columns
    gRNA = target.iloc[:, 0]
    gene = target.iloc[:, 1]
    LFC = target.iloc[:, 2]
    Q = target.iloc[:, 3]

    T_wt = pd.DataFrame({
        'gRNA': gRNA,
        'gene': gene,
        'LFC': LFC,
        'Z_zGE_LFC': Z_t,
        'Q': Q
    })

    return T_wt
