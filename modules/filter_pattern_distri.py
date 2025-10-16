import numpy as np
import pandas as pd 

def filter_pattern_distri(raw_ind, pat, n, st, en, step):
    """
    Parameters
    ----------
    raw_ind : pandas.DataFrame
        Input table containing gRNA information and read counts. 
        Must contain at least the following columns:
        ['gRNA', 'gene', 'T0_F', 'T0_M', 'T1_F', 'T1_M'].
    pat : str
        Pattern (substring) to filter gene names. Case-insensitive.
    n : int
        Number of characters from the start of the gene name to compare 
        against the pattern (e.g. 3 → 'chr' matches 'chr1', 'chrX', etc.).
    st : float
        Start of histogram range (e.g. -5).
    en : float
        End of histogram range (e.g. 5).
    step : float
        Bin width for histogram (e.g. 0.5).

    Returns
    -------
    indiv_Chr : pandas.DataFrame
        Subset of input rows where the gene name matches the pattern.
    indiv_noChr : pandas.DataFrame
        Subset of input rows that do not match the pattern.
    ind_out : list of int
        Row indices (in `raw_ind`) of matching entries.
    ind_in : list of int
        Row indices (in `raw_ind`) of non-matching entries.
    num_out_in : int
        Number of matching rows.
    bins : numpy.ndarray
        Array of histogram bin edges.
    hissFM : numpy.ndarray
        Raw histogram counts of summed LFCs (female + male).
    percFM : numpy.ndarray
        Normalized histogram values in percent (summing to 100%).
    LFC_1 : numpy.ndarray
        Log2 fold-change values for females: log2((T1_F + 1e-6) / (T0_F + 1e-6)).
    LFC_2 : numpy.ndarray
        Log2 fold-change values for males: log2((T1_M + 1e-6) / (T0_M + 1e-6))
    s_chr : numpy.ndarray
        Sum of counts across all four columns ['T0_F', 'T0_M', 'T1_F', 'T1_M'].
    st1 : float
        Copy of the input `st` value (for downstream consistency).
    en1 : float
        Copy of the input `en` value (for downstream consistency).
    T_lfc_pat : pandas.DataFrame
        Table containing gRNA ID, gene name, and LFC values per replicate:
        columns ['gRNA_chr', 'genes_chr', 'lfc1', 'lfc2'].
    """

    # 1. Extract genes and gRNAs
    genes = raw_ind.iloc[:, 1].astype(str)
    gRNA = raw_ind.iloc[:, 0]

    # 2. Apply pattern filter using only the first n characters
    def matches_pattern(g):
        return g[:n].lower() == pat.lower()
    
    matches = genes.apply(matches_pattern)
    ind_out = matches[matches].index  # pattern-matching rows
    ind_in = matches[~matches].index  # rest

    indiv_Chr = raw_ind.loc[ind_out]
    indiv_noChr = raw_ind.loc[ind_in]
    num_out_in = len(ind_out)

    # 3. Compute LFCs from columns 3 to 6
    # Assume: col3 = rep1_T1, col4 = rep1_T2, col5 = rep2_T1, col6 = rep2_T2
    # Correct LFCs: T1_F / T0_F and T1_M / T0_M
    counts = indiv_Chr[['T0_F', 'T0_M', 'T1_F', 'T1_M']].astype(float)

    LFC_1 = np.log2((counts['T1_F'] + 1e-6) / (counts['T0_F'] + 1e-6))
    LFC_2 = np.log2((counts['T1_M'] + 1e-6) / (counts['T0_M'] + 1e-6))

    LFC_sum = LFC_1 + LFC_2

    bins = np.arange(st, en + step, step)
    hissFM, _ = np.histogram(LFC_sum, bins)
    percFM = 100 * hissFM / hissFM.sum() if hissFM.sum() > 0 else np.zeros_like(hissFM)

    # Optional: total raw counts (for reference, not returned)
    s_chr = counts.sum().values
    st1, en1 = st, en

    # 4. Construct LFC table
    T_lfc_pat = pd.DataFrame({
        'gRNA_chr': gRNA.loc[ind_out].values,
        'genes_chr': genes.loc[ind_out].values,
        'lfc1': LFC_1.values,
        'lfc2': LFC_2.values
    })

    return (
        indiv_Chr, indiv_noChr,
        ind_out.tolist(), ind_in.tolist(), num_out_in,
        bins, hissFM, percFM,
        LFC_1.values, LFC_2.values,
        s_chr, st1, en1, T_lfc_pat
    )
