import pandas as pd
import numpy as np

def make_LFC_Z_MZ_tables_two(zGE_12: pd.DataFrame,
                             MZ12: np.ndarray,
                             Z12: np.ndarray):
    """
    Build per‑replicate tables that append Z‑score and MZ‑score columns
    to the zGE LFC data.

    Parameters
    ----------
    zGE_12 : pd.DataFrame
        Columns order: 0‑gRNA, 1‑gene, 2‑lfc1, 3‑lfc2
    MZ12   : ndarray, shape (n_rows, 2)
        MZ scores for replicate‑1 (col‑0) and replicate‑2 (col‑1)
    Z12    : ndarray, shape (n_rows, 2)
        Z  scores for replicate‑1 (col‑0) and replicate‑2 (col‑1)

    Returns
    -------
    T_zGE_1 : DataFrame  (replicate‑1)  columns: gRNA, gene, LFC, Z, MZ
    T_zGE_2 : DataFrame  (replicate‑2)  columns: gRNA, gene, LFC, Z, MZ
    """

    # shared identifiers
    gRNA = zGE_12.iloc[:, 0].values
    gene = zGE_12.iloc[:, 1].values

    # ---------- replicate‑1 ----------
    T_zGE_1 = pd.DataFrame({
        "gRNA": gRNA,
        "gene": gene,
        "LFC":  zGE_12.iloc[:, 2].astype(float).values,   # lfc1
        "Z":    Z12[:, 0].astype(float),
        "MZ":   MZ12[:, 0].astype(float),
    })

    # ---------- replicate‑2 ----------
    T_zGE_2 = pd.DataFrame({
        "gRNA": gRNA,
        "gene": gene,
        "LFC":  zGE_12.iloc[:, 3].astype(float).values,   # lfc2
        "Z":    Z12[:, 1].astype(float),
        "MZ":   MZ12[:, 1].astype(float),
    })

    return T_zGE_1, T_zGE_2

