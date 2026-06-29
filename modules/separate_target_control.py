from .filter_pattern_distri import filter_pattern_distri
from .zGE_target_distri import zGE_target_distri
from .compute_hiss_LFC_rep12 import compute_hiss_LFC_rep12
from .distri_target_contr_plots_all import distri_target_contr_plots_all
import numpy as np

def separate_target_control(d, st, en, step, T_norm_indiv, zGE, pat1, pat2, output_dir):
    """
    Parameters:
    - d: dataset used in zGE_target_distri
    - st, en, step: for histogram binning
    - T_norm_indiv: full normalized gRNA table
    - zGE: list (table-like) of zGE genes (first column = gene names)
    - pat1: pattern for intergenic control (e.g., 'chr')
    - pat2: pattern for NT control (e.g., 'Non')

    Returns:
    - T_target: target genes (excluding controls + zGE)
    - T_target_zGE: target + zGE genes
    - T_lfc_chr: intergenic control genes
    - T_lfc_nt: NT control genes
    - T_zGE: genes from zGE list (sorted)
    - bin: histogram bins
    - All histograms (his*, perc*) and gzn table
    """

    # Initialize outputs
    T_target = None
    T_target_zGE = None
    T_lfc_chr = None
    T_lfc_nt = None
    T_zGE = None
    gzn = None

    # -----------------------
    # 1. Separate intergenic (e.g., 'chr')
    pat = pat1
    n = len(pat)

    (
        indiv_Chr, indiv_noChr,
        ind_out1, ind_in1, num_out_in1,
        bin, hisFMi, percFMi,
        LFC_1i, LFC_2i,
        s_chr, st1, en1, T_lfc_chr
    ) = filter_pattern_distri(T_norm_indiv, pat, n, st, en, step)
    
    num_chr_rest = num_out_in1

    # -----------------------
    # 2. Separate NT (e.g., 'Non')
    pat = pat2
    n = len(pat)
    (
        indiv_NT, indiv_noChr_noNT,
        ind_out2, ind_in2, num_out_in2,
        bin, hisFMnt, percFMnt,
        LFC_1nt, LFC_2nt,
        s_nt, st1, en1, T_lfc_nt
    ) = filter_pattern_distri(indiv_noChr, pat, n, st, en, step)

    num_NT_rest = num_out_in2

    # -----------------------
    # 3. Separate zGE and target genes
    zGE_genes = zGE.iloc[:, 0].tolist()
    genes_nnt = indiv_noChr_noNT.iloc[:, 1].tolist()  # gene column

    (
        T_zGE, T_nzGE, bin,
        his1z, perc1z, his2z, perc2z, hisFMz, percFMz,
        his1t, perc1t, his2t, perc2t, hisFMt, percFMt,
        gzn
    ) = zGE_target_distri(d, genes_nnt, zGE_genes, indiv_noChr_noNT, st, en, step)

    # -----------------------
    # 4. Compute distributions for chr and NT controls
    bin, his1i,perc1i,his2i,perc2i,hisFMi,percFMi,LFC_1i,LFC_2i, st1i, en1i = compute_hiss_LFC_rep12(T_lfc_chr, st, en, step)
    bin, his1n,perc1n,his2n,perc2n,hisFMn,percFMn,LFC_1n,LFC_2n, st1n, en1n = compute_hiss_LFC_rep12(T_lfc_nt, st, en, step)

    # -----------------------
    # 5. Plot histograms (female and male)
    cond1 = "F"
    cond2 = "M"
    distri_target_contr_plots_all(bin, perc1z, perc1i, perc1t, perc1n, cond1, output_dir)
    distri_target_contr_plots_all(bin, perc2z, perc2i, perc2t, perc2n, cond2, output_dir)

    # -----------------------
    # 6. Final Outputs
    T_target = T_nzGE  # targets only
    T_target_zGE = indiv_noChr_noNT  # target + zGE
    T_zGE = T_zGE.sort_values(by="gene")  # sort alphabetically

    return (
        T_target, T_target_zGE, T_lfc_chr, T_lfc_nt, T_zGE,
        bin, his1z, perc1z, his2z, perc2z, hisFMz, percFMz,
        his1t, perc1t, his2t, perc2t, hisFMt, percFMt, gzn
    )