# Fowler_et_al

This pipeline is designed to find significant gRNA enrichment/depletion accounting for both FDR and FNR rates. This is achieved by calculating the likelihood of the observed gRNA fold change for each gene differing from the distribution of the set of true neutral control gRNAs in the screen.

## System Requirements

### Operating Systems
This package has been tested on:
- macOS: Tahoe (26.5.1)
- Linux: Ubuntu 22.04 LTS
- Windows: 11

### Software Dependencies
- Python 3.13.5

See required packages and versions in requirement.txt

### Hardware Requirements
Runs on a standard computer with enough RAM for in-memory operations (e.g., 8 GB). <br />
No non-standard hardware required.

## Installation

Clone the repository:
```
git clone https://github.com/PHJonesGroup/Fowler_etal_code.git
cd Fowler_et_al/
```

Create and activate a virtual environment:
```
python -m venv venv
source venv/bin/activate      # On Windows: venv\Scripts\activate
```

Install dependencies:
```
pip install -r requirements.txt
```

**Typical install time:** ~10 minutes on a standard desktop with a normal internet connection.

## Running the python script
To run on your own data:

1. Copy the example config: `cp config.yaml my_config.yaml`
2. Edit the parameters (see the table below).
3. Run: `python3 targeted_CRISPR.py --config my_config.yaml`

### Configuration parameters

| Parameter        | Description                                       | Example                  |
|------------------|---------------------------------------------------|--------------------------|
| `input_counts`   | Path to the counts CSV                            | `demo_data/WT_....csv`   |
| `input_controls` | Path to the control genes CSV                     | `demo_data/genes_....csv`|
| `output_dir`     | Directory where results are written              | `output`                 |
| `index_scheme`   | 2 = use both replicates; 1 = single replicate    | `2`                      |
| `d`              | Number of gRNAs per gene                          | `4`                      |
| `alf`            | Significance level (tail of the distribution)     | `0.06`                   |
| `st`            | Start of x axis     | `-10`                   |
| `end`            | End of x axis     | `10`                   |
| `step`            | Bin sizes     | `0.05`                   |
| `thr_lfchz`            | LFC threshold for volcano plot     | `1.7`                   |
| `pat1`, `pat2`   | String patterns used to filter genes              | `'chr'`, `'Non'`         |
| `cond1`, `cond2` | Condition labels              | `'F'`, `'M'`             |

### Pipeline Outputs
| File Name        | Description                                       |
|------------------|---------------------------------------------------|
| `gene_distribution_gRNA.png`   | Path to the counts CSV                            | 
| `gRNA_counts_normalisation.png`   | Path to the counts CSV                            | 
| `normalised_counts_mean_med.png`   | Path to the counts CSV                            | 
| `mean_med_T0_T1.png`   | Path to the counts CSV                            | 
| `distri_separate_target_controls_{condition_name1/2}_WTval.png`   | Path to the counts CSV                            | 
| `distri_target_controls_{condition_name1/2}_WTval.png`   | Path to the counts CSV                            | 
| `p_distri_targ_cont_{condition_name1/2}.png`   | Path to the counts CSV                            | 
| `distri_LFC_zGE_{condition_name1/2}.png`   | Path to the counts CSV                            | 
| `p_controls_rep1_rep2_zGE_{condition_name1+2}.png`   | Path to the counts CSV                            | 
| `target_distri_LFC_Z_corr_Target_genes.png`   | Path to the counts CSV                            | 
| `target_distri_LFC_Z_corr_Interg.png`   | Path to the counts CSV                            | 
| `target_distri_LFC_Z_corr_NT.png`   | Path to the counts CSV                            | 
| `volcano_gRNA_{condition_name1/2}.png`   | Path to the counts CSV                            | 
| `volcano_plot_interactive_{condition_name1/2}.html`   | Path to the counts CSV                            | 
| `T_NT_{condition_name1/2}_wt_gRNA.csv`   | Path to the counts CSV                            | 
| `T_intergenic_{condition_name1/2}_wt_gRNA.csv `   | Path to the counts CSV                            | 
| `T_intergenic_{condition_name1/2}_wt_gRNA.csv  `   | Path to the counts CSV                            | 
| `T_norm_T0T1_indiv_WT.csv`   | Path to the counts CSV                            | 
| `T_target_{condition_name}_wt_gRNA.csv`   | Path to the counts CSV                            | 


## Pipeline Steps
### Module 1
Normalise raw count data as proportion per million <br />
Counts number of gRNA per gene/ intergenic <br />
Compute mean and median of normalised counts <br />

### Module 2
Compute and analyse LFC distributions of targets and controls (target, intergenic and zero expressed genes) <br />
Compute q-values, critical values, mean and standard deviation for zero expressed genes controls. <br />
Compute z for targets and all controls, based on mean and standard deviation of zero expressed genes and controls. <br />
Determine significant gRNA enrichment/depletion accounting for both FDR and FNR rates <br />
