import pandas as pd

def norm_table_individ_WT_valid(norm_dat, raw_ind):
    """
    Normalize and create tables for WT validation for T0 and T1 timepoints.

    Inputs:
        norm_dat: numpy array or DataFrame with normalized data, expected shape (N,4)
                  Columns: [F_T0, F_T1, M_T0, M_T1]
        raw_ind: pandas DataFrame with columns:
                 ['sgRNA_name', 'gene', ..., others]

    Outputs:
        tab_norm_T0: pandas DataFrame with columns ['sgRNA_name', 'gene', 'T0_F', 'T0_M']
        tab_norm_T1: pandas DataFrame with columns ['sgRNA_name', 'gene', 'T1_F', 'T1_M']
    """
    # Check that norm_dat has 4 columns as expected
    if norm_dat.shape[1] == 4:
        # Extract T0 and T1 data 
        T0_dat = norm_dat[:, [0, 2]]  
        T1_WT_dat = norm_dat[:, [1, 3]]  
    else:
        raise ValueError("wrong format of data: expected 4 columns in norm_dat")

    # Extract sgRNA and gene columns from raw_ind DataFrame
    gRNA = raw_ind.iloc[:, 0]
    gene = raw_ind.iloc[:, 1]

    # Create pandas DataFrames for T0 and T1 normalized data
    tab_norm_T0 = pd.DataFrame({
        'sgRNA_name': gRNA,
        'gene': gene,
        'T0_F': T0_dat[:, 0],
        'T0_M': T0_dat[:, 1]
    })

    tab_norm_T1 = pd.DataFrame({
        'sgRNA_name': gRNA,
        'gene': gene,
        'T1_F': T1_WT_dat[:, 0],
        'T1_M': T1_WT_dat[:, 1]
    })

    return tab_norm_T0, tab_norm_T1