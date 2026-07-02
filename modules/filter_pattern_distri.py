import numpy as np
import pandas as pd

def filter_pattern_distri(raw_ind, pat, n, st, en, step, cond1, cond2):
    # 1. Extract genes and gRNAs
    genes = raw_ind.iloc[:, 1].astype(str)
    gRNA  = raw_ind.iloc[:, 0]

    # 2. Pattern filter on first n characters
    matches = genes.apply(lambda g: g[:n].lower() == pat.lower())
    ind_out = matches[matches].index
    ind_in  = matches[~matches].index

    indiv_Chr   = raw_ind.loc[ind_out]
    indiv_noChr = raw_ind.loc[ind_in]
    num_out_in  = len(ind_out)

    # 3. Find all replicate columns for each condition
    cond1_cols = [c for c in raw_ind.columns if c.startswith(f"{cond1}_")]
    cond2_cols = [c for c in raw_ind.columns if c.startswith(f"{cond2}_")]
    if not cond1_cols or not cond2_cols:
        raise ValueError(f"Missing columns for {cond1}_ or {cond2}_")

    # 4. One LFC per gRNA: log2( mean(cond1 reps) / mean(cond2 reps) )
    eps = 1e-6
    if len(indiv_Chr):
        c1_mean = indiv_Chr[cond1_cols].astype(float).mean(axis=1)
        c2_mean = indiv_Chr[cond2_cols].astype(float).mean(axis=1)
        LFC = np.log2((c1_mean + eps) / (c2_mean + eps)).values   # 1-D, one per gRNA
    else:
        LFC = np.empty(0)

    bins = np.arange(st, en + step, step)
    hissFM, _ = np.histogram(LFC, bins)
    percFM = 100 * hissFM / hissFM.sum() if hissFM.sum() > 0 else np.zeros_like(hissFM, dtype=float)

    all_cols = cond1_cols + cond2_cols
    s_chr = indiv_Chr[all_cols].astype(float).sum().values
    st1, en1 = st, en

    # 5. LFC table: ids + single lfc column
    T_lfc_pat = pd.DataFrame({
        'gRNA':  gRNA.loc[ind_out].values,
        'genes': genes.loc[ind_out].values,
        'lfc':   LFC,
    })

    return (
        indiv_Chr, indiv_noChr,
        ind_out.tolist(), ind_in.tolist(), num_out_in,
        bins, hissFM, percFM,
        LFC,
        s_chr, st1, en1, T_lfc_pat
    )