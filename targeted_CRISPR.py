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
ssc = config["ssc"]
alf = config["alf"]
st = config["st"]
en = config["en"]
step = config["step"]
start_end = config["start_end"]
pat1 = config["pat1"]
pat2 = config["pat2"]
n = config["n"]
cond1 = config["cond1"]
cond2 = config["cond2"]
condz = config["condz"]
indiv = config["indiv"]
cond_id = config["cond_id"]
cond11 = config["cond11"]
cond12 = config["cond12"]
condz1 = config["condz1"]
condz2 = config["condz2"]
sfdr_corr = config["sfdr_corr"]
thr_lfchz = config["thr_lfchz"]
thr_lfcdz = config["thr_lfcdz"]

# ============================ MODULE 1 ============================
raw_ind = pd.read_csv(input_counts)
si_input = raw_ind.shape

T_vert = raw_ind
gRNA = raw_ind.iloc[:, 0].values  # first column (gRNA names)
gene = raw_ind.iloc[:, 1].values  # second column (gene names)

# -------------------- 1.1 --------------------
# Counts how many times each gene name/intergenic name occurs in the raw count file & gives number of gRNA per gene/ intergenic
gene_names = gene
ggenes, gg, ind_gn = count_vertical_names.count_vertical_names(gene_names)
plot_grna_distribution.plot_grna_distribution(gg,d, output_dir)

# -------------------- 1.2 --------------------
# Removes two- and three- gRNA block gene
d = 4 if max(gg) == 4 else max(gg)

if d == 4:
    (Genes_1, Genes_2, Genes_3, Genes_4,
     gn_1, gn_2, gn_3, gn_4,
     wt_1, wt_2, wt_3, wt_4,
     d, nums,
     ind4, ind3, ind2, ind1) = separate_fours_threes_twos_ones_genes_gRNAs.separate_fours_threes_twos_ones_genes_gRNAs(T_vert, ggenes, gene_names, ind_gn, gg)
    num_gRNA = nums

if nums[1] > 0:
    Genes_3_head = Genes_3[:3]
else:
    print(Genes_3)

if nums[2] > 0:
    Genes_2_head = Genes_2[:3]
else:
    print(Genes_2)

raw_ind_4 = T_vert.iloc[ind4, :]
raw_ind_1 = T_vert.iloc[ind1, :]

raw_ind_f = pd.concat([raw_ind_4, raw_ind_1], ignore_index=True)

si_4_1 = raw_ind_f.shape

# -------------------- 1.3 --------------------
# Normalise (FPKM) counts as percentage within a column
norm_dat = normalise_prop.normalise_prop(ssc, raw_ind_f, output_dir) 

# -------------------- 1.4 --------------------
# Format normalised counts, compute mean and median of normalised counts
tab_norm_T0, tab_norm_T1 = norm_table_individ_WT_valid.norm_table_individ_WT_valid(norm_dat, raw_ind_f)
me_med_T0_T1 = plot_me_med.plot_me_med(tab_norm_T0, tab_norm_T1, output_dir)

T_norm_WT = pd.concat([tab_norm_T0, tab_norm_T1.iloc[:, 2:4]], axis=1)
T_norm_WT.to_csv(os.path.join(output_dir, 'T_norm_T0T1_indiv_WT.csv'), index=False)

norm_tab = T_norm_WT.copy()
column_labels = norm_tab.columns.tolist()
num_gRNA = d

T_norm_indiv, norm_dat, me_med_nd = format_data.format_data(index_scheme, indiv, norm_tab, output_dir)

# ============================ MODULE 2 ============================

# Load zGE dataset
zGE = pd.read_csv(input_controls)
zGE = zGE.drop_duplicates() 

# Check full duplicate rows
dup_rows = T_norm_indiv[T_norm_indiv.duplicated()]
if not dup_rows.empty:
    print(dup_rows.head())

# Assuming gRNA is in the first column
gRNA_col = T_norm_indiv.columns[0]
dup_gRNAs = T_norm_indiv[gRNA_col][T_norm_indiv[gRNA_col].duplicated()]
if not dup_gRNAs.empty:
    print(dup_gRNAs.value_counts().head())

# -------------------- 2.1 --------------------
# Separate and show controls, NTs, zGE, and normalised counts
# Distribution of:
# (1) Controls: intergenic genes
# (2) NT = non-targeted genes
# (3) zGE = zero expressed genes
# (4) Normalised targeted genes

(
    T_target, T_target_zGE_counts, T_lfc_chr, T_lfc_nt, T_zGE,
    bin, his1z, perc1z, his2z, perc2z, hisFMz, percFMz,
    his1t, perc1t, his2t, perc2t, hisFMt, percFMt, num_gzt
) = separate_target_control.separate_target_control(
    d=d, st=st, en=en, step=step,
    T_norm_indiv=T_norm_indiv, zGE=zGE,
    pat1=pat1, pat2=pat2, n=n, output_dir=output_dir
)

# -------------------- 2.2 --------------------
# Choose zGE controls, compute Z/MZ/crit‑LR, implement q, adjust controls for Z

his12t = np.column_stack([his1t, his2t])

(
    crit_LR12,
    me_sd12,
    med_mad12,
    binn,
    p_cont12,
    hiss_cont12,
    p_targ12,
    T_zGE_1,
    T_zGE_2,
) = CTR_stats_zGE.CTR_stats_zGE(
    alf,
    st,
    en,
    step,
    T_zGE,          # DataFrame with gRNA / gene / lfc1 / lfc2
    his12t,         
    condz1,
    condz2,
    cond1,          
    cond2,          
    condz,
    output_dir
)

# -------------------- 2.3 --------------------
# Get recalibrated p/q values for gRNAs
T_vert_q = p_control_target_implement.p_control_target_implement(T_target, T_zGE, bin, p_cont12)

# -------------------- 2.4 --------------------
# Compute Z LFC for two replicas rep1 rep2, for targets only

LFC_t1 = T_vert_q.iloc[:, 2].to_numpy()
LFC_t2 = T_vert_q.iloc[:, 3].to_numpy()

Z_t1, Z_t2, bin, perc_t1, perc_zt1, perc_t2, perc_zt2 = computeZ.computeZ(
    st, en, step, LFC_t1, LFC_t2, me_sd12
)

plot_histograms.plot_histograms(bin, perc_t1, perc_zt1, perc_t2, perc_zt2, "Target_genes", output_dir)

# make tables gRNA-based
T_t1, T_t2 = make_tables_Z_two.make_tables_Z_two(T_vert_q, Z_t1, Z_t2, cond1, cond2)

# Module for Z any set: intergenic, NT, etc
# Intergenic
T_any = T_lfc_chr.copy()
T_z_q_chr1, T_z_q_chr2, binn_chr, perc_t1_chr, perc_zt1_chr, perc_t2_chr, perc_zt2_chr = z_p_CTR_any.z_p_CTR_any(
    st, en, step, T_any, binn, p_cont12, me_sd12, cond1, cond2
)

plot_histograms.plot_histograms(binn_chr, perc_t1_chr, perc_zt1_chr, perc_t2_chr, perc_zt2_chr, cond11, output_dir)

# NT control
T_any = T_lfc_nt.copy()
T_z_q_nt1, T_z_q_nt2, binn_nt, perc_t1_nt, perc_zt1_nt, perc_t2_nt, perc_zt2_nt = z_p_CTR_any.z_p_CTR_any(
    st, en, step, T_any, binn, p_cont12, me_sd12, cond1, cond2
)

plot_histograms.plot_histograms(binn_nt, perc_t1_nt, perc_zt1_nt, perc_t2_nt, perc_zt2_nt, cond12, output_dir)

# -------------------- 2.5 --------------------
# Find extreme sets (enriched/ depleted) and volcano for gRNA and gene

# rep1
thr_lfch = crit_LR12[0, 1]  
thr_lfcd = crit_LR12[0, 0]

thrLFC_d_h = [thr_lfcd, thr_lfch]


(T_gRNA_1, T_lfc_z_q_med_1, T_lfc_z_q_me_1, T_LFC_1, T_Q_1, T_Z_1, 
 indha_1, indda_1, indhaz_1, inddaz_1) = volcano_gRNA_gene_hits_wt.volcano_gRNA_gene_hits_wt(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_t1, cond1, output_dir)

(T_gRNA_1, T_lfc_z_q_med_1, T_lfc_z_q_me_1, T_LFC_1, T_Q_1, T_Z_1, 
 indha_1, indda_1, indhaz_1, inddaz_1) = volcano_gRNA_gene_hits_wt_interactive.volcano_gRNA_gene_hits_wt_interactive(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_t1, cond1, output_dir)


num_hd_LFC_Z_1 = [len(indha_1), len(indda_1), len(indhaz_1), len(inddaz_1)]

# rep2
numeric_cols = ['LFC', 'Z_zGE_LFC', 'Q']
T_gene = T_t2.groupby('gene')[numeric_cols].median().reset_index()

gene_counts = T_t2['gene'].value_counts()

thr_lfch = crit_LR12[1, 1]  
thr_lfcd = crit_LR12[1, 0]

(T_gRNA_2, T_lfc_z_q_med_2, T_lfc_z_q_me_2, T_LFC_2, T_Q_2, T_Z_2,
 indha_2, indda_2, indhaz_2, inddaz_2) = volcano_gRNA_gene_hits_wt.volcano_gRNA_gene_hits_wt(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_t2, cond2, output_dir)

(T_gRNA_2, T_lfc_z_q_med_2, T_lfc_z_q_me_2, T_LFC_2, T_Q_2, T_Z_2,
 indha_2, indda_2, indhaz_2, inddaz_2) = (
    volcano_gRNA_gene_hits_wt_interactive.volcano_gRNA_gene_hits_wt_interactive(
        alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_t2, cond2, output_dir
    )
)

num_hd_LFC_Z_2 = [len(indha_2), len(indda_2), len(indhaz_2), len(inddaz_2)]

# Saving results
T_z_q_chr1.to_csv(os.path.join(output_dir, 'T_intergenic_F_wt_gRNA.csv'), index=False)
T_z_q_chr2.to_csv(os.path.join(output_dir, 'T_intergenic_M_wt_gRNA.csv'), index=False)
T_z_q_nt1.to_csv(os.path.join(output_dir, 'T_NT_F_wt_gRNA.csv'), index=False)
T_z_q_nt2.to_csv(os.path.join(output_dir, 'T_NT_M_wt_gRNA.csv'), index=False)
T_t2.to_csv(os.path.join(output_dir, 'T_targetM_wt_gRNA.csv'), index=False)
