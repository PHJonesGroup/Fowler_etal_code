import pandas as pd
import numpy as np

def make_LFC_Z_MZ_tables_two(zGE: pd.DataFrame,
                             MZ: np.ndarray,
                             Z: np.ndarray):
    """
    Build table that append Z‑score and MZ‑score columns
    to the zGE LFC data.

    Parameters
    ----------
    zGE : pd.DataFrame
        Columns order: 0‑gRNA, 1‑gene, 2‑lfc
    MZ   : ndarray, shape (n_rows, 2)
        MZ scores for replicate‑1 (col‑0)
    Z    : ndarray, shape (n_rows, 2)
        Z  scores for replicate‑1 (col‑0) 

    Returns
    -------
    T_zGE : DataFrame columns: gRNA, gene, LFC, Z, MZ
    """

    # shared identifiers
    gRNA = zGE.iloc[:, 0].values
    gene = zGE.iloc[:, 1].values

    # ---------- replicate‑1 ----------
    T_zGE = pd.DataFrame({
        "gRNA": gRNA,
        "gene": gene,
        "LFC":  zGE.iloc[:, 2].astype(float).values,
        "Z":    np.asarray(Z,  dtype=float),
        "MZ":   np.asarray(MZ, dtype=float),
    })

    return T_zGE

