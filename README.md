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
- contourpy==1.3.2
- cycler==0.12.1
- fonttools==4.59.0
- Jinja2==3.1.6
- kiwisolver==1.4.8
- MarkupSafe==3.0.2
- matplotlib==3.10.3
- narwhals==1.48.1
- numpy==2.3.2
- packaging==25.0
- pandas==2.3.1
- pillow==11.3.0
- plotly==6.2.0
- pyan3==1.2.0
- pydeps==3.0.1
- pyparsing==3.2.3
- python-dateutil==2.9.0.post0
- pytz==2025.2
- PyYAML==6.0.2
- scipy==1.16.0
- seaborn==0.13.2
- six==1.17.0
- snakeviz==2.2.2
- stdlib-list==0.11.1
- tornado==6.5.1
- tzdata==2025.2

### Hardware Requirements
Runs on a standard computer with enough RAM for in-memory operations (e.g., 8 GB). <br />
No non-standard hardware required.

## Installation

### Install

Clone the repository:
\```bash
git clone https://github.com/PHJonesGroup/Fowler_etal_code.git
cd Fowler_et_al/
\```

Create and activate a virtual environment:
\```bash
python -m venv venv
source venv/bin/activate      # On Windows: venv\Scripts\activate
\```

Install dependencies:
\```bash
pip install -r requirements.txt
\```

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
Before running the script, ensure that the config.yaml is populated.
```
python3 targeted_CRISPR.py   
```
