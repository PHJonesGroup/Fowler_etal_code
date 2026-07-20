# Fowler_et_al

This pipeline is designed to find significant gRNA enrichment/depletion accounting for both FDR and FNR rates. This is achieved by calculating the likelihood of the observed gRNA fold change for each gene differing from the distribution of the set of true neutral control gRNAs in the screen.

## System Requirements

### Operating Systems
This package has been tested on:
- macOS: Tahoe (26.5.1)
- Linux: Ubuntu 22.04 LTS

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
cd Fowler_etal_code/
```

Create and activate a virtual environment:
```
python -m venv myenv
source myenv/bin/activate
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
| `input_controls` | Path to the control genes CSV (optional)                  | `demo_data/genes_....csv`|
| `output_dir`     | Directory where results are written              | `output`                 |
| `rep`   | Number of replicates per condition    | `4`                      |
| `keep_counts`   | Number of gRNAs expected per genes    | `[1, 4, 10]`                      |
| `control`   |  Baseline against which targeted genes distribution is compared against. Choose between Intergenic, Non-targetting or zGE   | `Non-targetting`                      |
| `alf`            | Significance level (tail of the distribution)     | `0.06`                   |
| `pat1`, `pat2`   | String patterns used to filter genes              | `'chr'`, `'Non'`         |
| `cond1`, `cond2` | Condition labels              | `'T1'`, `'T0'`             |
| `st`            | Start of x axis     | `-10`                   |
| `end`            | End of x axis     | `10`                   |
| `step`            | Bin sizes     | `0.05`                   |
| `thr_lfchz`            | LFC threshold for volcano plot     | `1.7`     |

### Outputs
| File Name        | Description                                       |
|------------------|---------------------------------------------------|
| `distri_separate_target_controls_{condition_name1_2}.png`   | Individual LFC distributions of non-targetting, zero expressed, intergenic and target genes                        | 
| `distri_target_controls_{condition_name1_2}.png`   | Overlay of LFC distributions of non-targetting, zero expressed, intergenic and target genes                            | 
| `distri_zGE_{condition_name1_2}.png`   | Zero Gene Expression distribution of LFC, Z-corrected LFC and MZ-corrected LFC (optional if zGE list given)                          | 
| `gene_distribution_gRNA.png`   | Number of gRNAs attributed to each gene                      | 
| `gRNA_counts_normalisation.png`   | Total raw (top) and normalised (bottom) read count summed over all gRNAs in each condition                            | 
| `normalised_counts_mean_med.png`   | Mean versus median of the normalised gRNA counts for each of condition                           | 
| `p_controls_zGE_{condition_name1_2}.png`   | P-value curve for zero expressed gRNAs (optional if zGE list given)                        | 
| `p_distri_targ_cont_{condition_name1_2}.png`   | Control-calibrated p-value curves of target vs control gRNAs                           |  
| `normalised_counts.csv`   | Normalised counts in csv format                            | 

## Pipeline Steps
### Module 1
Normalise raw count data as proportion per million <br />
Counts number of gRNA per gene/ intergenic <br />
Compute mean and median of normalised counts <br />

### Module 2
Compute and analyse LFC distributions of targets and controls (target, intergenic and zero expressed genes) <br />
Compute q-values, critical values, mean and standard deviation for zero expressed genes controls. <br />
Compute z for targets and all controls, based on mean and standard deviation of zero expressed genes and controls. <br />
