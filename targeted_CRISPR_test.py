import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import yaml
from pathlib import Path
import argparse

from modules import count_vertical_names, plot_grna_distribution
from modules import separate_fours_threes_twos_ones_genes_gRNAs
from modules import normalise_prop
from modules import norm_table_individ_WT_valid, plot_me_med
from modules import find_columns_indiv, format_data

from modules import compute_hiss_LFC_rep12, distri_target_contr_plots_all, filter_pattern_distri, make_histo_LFC, separate_target_control, zGE_target_distri
from modules import compute_p_critLFC, CTR_stats_zGE, make_histo_crit_stats, make_histo_vec_rep12, make_LFC_Z_MZ_tables_two, med_mad_MZNP_2, q_val_frequentist_critical
from modules import implement_p_control_indiv_gRNA, p_control_any_implement, p_control_target_implement, plot_histograms
from modules import computeZ, make_tables_Z_two, z_p_CTR_any
from modules import perGene_4_hits_med_horiz, perGene_4_med_horiz, separate_fours_threes_twos_genes, volcano_gRNA_gene_hits_wt_interactive, volcano_gRNA_gene_hits_wt
from modules import separate_genes_by_grna_count

def parse_args():
    parser = argparse.ArgumentParser(description="Run the pipeline.")
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to the YAML config file (default: config.yaml)",
    )
    return parser.parse_args()

def load_config(path):
    if not os.path.exists(path):
        sys.exit(f"Config error: file not found: {path}")
    with open(path, "r") as f:
        try:
            return yaml.safe_load(f)
        except yaml.YAMLError as e:
            sys.exit(f"Config error: could not parse {path}\n{e}")

args = parse_args()
config = load_config(args.config)

input_counts = Path(config["input_counts"])
input_controls = Path(config["input_controls"])
output_dir = Path(config["output_dir"])
index_scheme = config["index_scheme"]
d = config["d"]
rep = config["rep"]
keep_counts = config["keep_counts"]
control = config["control"]
alf = config["alf"]
st = config["st"]
en = config["en"]
step = config["step"]
pat1 = config["pat1"]
pat2 = config["pat2"]
cond1 = config["cond1"]
cond2 = config["cond2"]
thr_lfchz = config["thr_lfchz"]

condz   = f'zGE_{cond1}{cond2}'
cond11  = 'Intergenic'
cond12  = 'Non_Targetting'
indiv   = f'{cond1}{cond2}'  
cond_id = f'{cond1}_{cond2}'
thr_lfcdz = -thr_lfchz
sfdr_corr = 0.0001 # constant to avoid log10(0)
ssc     = 1 # small sample correction count


# ============================ MODULE 1 ============================
raw_ind = pd.read_csv(input_counts)
si_input = raw_ind.shape

T_vert = raw_ind
gRNA = raw_ind.iloc[:, 0].values  # first column (gRNA names)
gene = raw_ind.iloc[:, 1].values  # second column (gene names)

# filter T_vert for columns of interest (cond1 and cond2)
keep = list(T_vert.columns[:2]) + [
    c for c in T_vert.columns[2:]
    if c.startswith(f"{cond1}_") or c.startswith(f"{cond2}_")
]
T_vert = T_vert[keep]

# -------------------- 1.1 --------------------
# Counts how many times each gene name/intergenic name occurs in the raw count file & gives number of gRNA per gene/ intergenic
gene_names = gene
ggenes, gg, ind_gn = count_vertical_names.count_vertical_names(gene_names)
plot_grna_distribution.plot_grna_distribution(gg,d, output_dir)

# -------------------- 1.2 --------------------
# Removes genes with incorrect number of gRNAs
groups, nums = separate_genes_by_grna_count.separate_genes_by_grna_count(T_vert, ggenes, gene_names, ind_gn, gg, rep)
ind_keep = sum((groups[c]["ind"] for c in keep_counts if c in groups), [])
raw_ind_f = T_vert.iloc[ind_keep, :].reset_index(drop=True)

# -------------------- 1.3 --------------------
# Normalise (FPKM) counts as percentage within a column
norm_dat = normalise_prop.normalise_prop(ssc, raw_ind_f, rep, output_dir) 
# -------------------- 1.4 --------------------
# Format normalised counts, compute mean and median of normalised counts
tab_norm_T0, tab_norm_T1 = norm_table_individ_WT_valid.norm_table_individ_WT_valid(norm_dat, raw_ind_f,  cond1, cond2, rep)
me_med_T0_T1 = plot_me_med.plot_me_med(tab_norm_T0, tab_norm_T1, output_dir)

T_norm_WT = pd.concat([tab_norm_T0, tab_norm_T1.iloc[:, 2:(2 + 2*rep)]], axis=1)
T_norm_WT.to_csv(os.path.join(output_dir, 'normalised_counts.csv'), index=False)

norm_tab = T_norm_WT.copy()
column_labels = norm_tab.columns.tolist()
num_gRNA = d

# ============================ MODULE 2 ============================
T_norm_indiv = T_norm_WT.copy()

# Load zGE dataset
zGE = pd.read_csv(input_controls)
zGE = zGE.drop_duplicates() 

# Check full duplicate rows
dup_rows = T_norm_indiv[T_norm_indiv.duplicated()]

# Assuming gRNA is in the first column
gRNA_col = T_norm_indiv.columns[0]
dup_gRNAs = T_norm_indiv[gRNA_col][T_norm_indiv[gRNA_col].duplicated()]

# -------------------- 2.1 --------------------
# Separate and show controls, NTs, zGE, and normalised counts
# Distribution of:
# (1) Controls: intergenic genes
# (2) NT = non-targeted genes
# (3) zGE = zero expressed genes
# (4) Normalised targeted genes

(
    T_target, T_target_zGE, T_lfc_chr, T_lfc_nt, T_zGE,
        bin, hisz, percz, hisFMz, percFMz,
        hist, perct, hisFMt, percFMt, gzn
) = separate_target_control.separate_target_control(
    d, st, en, step,
    T_norm_indiv, zGE,
    pat1, pat2, cond1, cond2, output_dir
)
print(T_target)
# -------------------- 2.2 --------------------
# Choose zGE controls, compute Z/MZ/crit‑LR, implement q, adjust controls for Z
if control == "Intergenic":
    (
        crit_LR, 
        me_sd, 
        med_mad, 
        binn, 
        p_cont, 
        hiss_cont, 
        p_targ, 
        T_zGE
    ) = CTR_stats_zGE.CTR_stats_zGE(
        alf,
        st,
        en,
        step,
        T_lfc_chr,          # DataFrame with gRNA / gene / lfc
        hist,         
        cond1,          
        cond2,          
        condz,
        control,
        output_dir
    )
elif control == "zGE":
    (
        crit_LR, 
        me_sd, 
        med_mad, 
        binn, 
        p_cont, 
        hiss_cont, 
        p_targ, 
        T_zGE
    ) = CTR_stats_zGE.CTR_stats_zGE(
        alf,
        st,
        en,
        step,
        T_zGE,          # DataFrame with gRNA / gene / lfc
        hist,         
        cond1,          
        cond2,          
        condz,
        output_dir
    )
else:
    print("Invalid input for control given.")

# -------------------- 2.3 --------------------
# Get recalibrated p/q values for gRNAs
T_vert_q = p_control_target_implement.p_control_target_implement(T_target, T_zGE, bin, p_cont)
# -------------------- 2.4 --------------------
# Compute Z LFC for two replicas rep1 rep2, for targets only

LFC_t = T_vert_q.iloc[:, 2].to_numpy()
Z_t, bin, perc_t, perc_zt,= computeZ.computeZ(
    st, en, step, LFC_t, me_sd
)
plot_histograms.plot_histograms(bin, perc_t, perc_zt, cond1, cond2, "Target_genes", output_dir)
print(T_vert_q)

# make tables gRNA-based
T_t = make_tables_Z_two.make_tables_Z_two(T_vert_q, Z_t)

# Module for Z any set: intergenic, NT, etc
if control == "zGE":
    # Intergenic
    T_any = T_lfc_chr.copy()
    T_z_q_chr1, binn_chr, perc_t1_chr, perc_zt1_chr = z_p_CTR_any.z_p_CTR_any(
        st, en, step, T_any, binn, p_cont, me_sd, cond1, cond2
    )

    plot_histograms.plot_histograms(binn_chr, perc_t1_chr, perc_zt1_chr, cond1, cond2, cond11, output_dir)

# NT control
T_any = T_lfc_nt.copy()
T_z_q_nt1, binn_nt, perc_t1_nt, perc_zt1_nt = z_p_CTR_any.z_p_CTR_any(
    st, en, step, T_any, binn, p_cont, me_sd, cond1, cond2
)

plot_histograms.plot_histograms(binn_nt, perc_t1_nt, perc_zt1_nt, cond1, cond2, cond12, output_dir)

# -------------------- 2.5 --------------------
# Find extreme sets (enriched/ depleted) and volcano for gRNA and gene

thr_lfch = crit_LR[1]  
thr_lfcd = crit_LR[0]

thrLFC_d_h = [thr_lfcd, thr_lfch]


(T_gRNA_1, T_lfc_z_q_med_1, T_lfc_z_q_me_1, T_LFC_1, T_Q_1, T_Z_1, 
 indha_1, indda_1, indhaz_1, inddaz_1) = volcano_gRNA_gene_hits_wt.volcano_gRNA_gene_hits_wt(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_t, cond1, output_dir)

(T_gRNA_1, T_lfc_z_q_med_1, T_lfc_z_q_me_1, T_LFC_1, T_Q_1, T_Z_1, 
 indha_1, indda_1, indhaz_1, inddaz_1) = volcano_gRNA_gene_hits_wt_interactive.volcano_gRNA_gene_hits_wt_interactive(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_t, cond1, output_dir)


num_hd_LFC_Z_1 = [len(indha_1), len(indda_1), len(indhaz_1), len(inddaz_1)]



# Saving results
if control == "zGE":
    T_z_q_chr1.to_csv(os.path.join(output_dir, 'intergenic_gRNA.csv'), index=False)

T_z_q_nt1.to_csv(os.path.join(output_dir, 'NT_gRNA.csv'), index=False)
T_t.to_csv(os.path.join(output_dir, 'target_gRNA.csv'), index=False)
