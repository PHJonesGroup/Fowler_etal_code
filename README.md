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

## Inputs
raw counts <br />
list of controls/ zero expressed genes (zGE) names <br />

## Pipeline Steps
### Module 1
Normalise raw count data as proportion per million <br />
Counts number of gRNA per gene/ intergenic <br />
Compute mean and median of normalised counts <br />

### Module 2
Compute and analyse LFC distributions of Targets and Controls (target, intergenic and zGE) <br />
Compute q-values, critical values, mean and std for zGE controls. <br />
Compute z for Targets and all CTRL, based on mean and std of zGE CTRL <br />
Determine significant gRNA enrichment/depletion accounting for both FDR and FNR rates <br />

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
| `pat1`, `pat2`   | String patterns used to filter genes              | `'chr'`, `'Non'`         |
| `cond1`, `cond2` | Condition labels (e.g. F / M)                     | `'F'`, `'M'`             |
