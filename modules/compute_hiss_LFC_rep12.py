import numpy as np
import pandas as pd
from .make_histo_LFC import make_histo_LFC

def compute_hiss_LFC_rep12(chr_lfc_12, st, en, step):
    """
    Compute histograms (distributions) for individual gRNAs, for both replicates and the average.
    
    Parameters:
        chr_lfc_12 (pd.DataFrame): DataFrame with columns:
            - gRNA_chr (str)
            - genes_chr (str)
            - lfc1 (float): replicate 1
            - lfc2 (float): replicate 2
        st (float): start of histogram range
        en (float): end of histogram range
        step (float): bin width

    Returns:
        binn: histogram bin edges
        hiss1: histogram counts for replicate 1
        perc1: histogram percentages for replicate 1
        hiss2: histogram counts for replicate 2
        perc2: histogram percentages for replicate 2
        hissFM: histogram counts for average
        percFM: histogram percentages for average
        LFC_1: replicate 1 LFCs (numpy array)
        LFC_2: replicate 2 LFCs (numpy array)
        st1: min(LFCs) - 1
        en1: max(LFCs) + 1
    """
    
    # Extract values
    LFC_1 = chr_lfc_12.iloc[:, 2].astype(float).values
    LFC_2 = chr_lfc_12.iloc[:, 3].astype(float).values

    # Compute new bounds
    st1 = min(LFC_1.min(), LFC_2.min()) - 1
    en1 = max(LFC_1.max(), LFC_2.max()) + 1

    # Compute average LFC (mean of both replicates)
    meLFC = np.mean(np.vstack((LFC_1, LFC_2)), axis=0)

    # Histograms
    binn, hissFM, percFM = make_histo_LFC(step, meLFC, st, en)
    _, hiss1, perc1 = make_histo_LFC(step, LFC_1, st, en)
    _, hiss2, perc2 = make_histo_LFC(step, LFC_2, st, en)

    return binn, hiss1, perc1, hiss2, perc2, hissFM, percFM, LFC_1, LFC_2, st1, en1

