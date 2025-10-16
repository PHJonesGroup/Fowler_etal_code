# Fowler_et_al

 This pipeline is designed to find significant gRNA enrichment/depletion accounting for both FDR and FNR rates. This is achieved by calculating the likelihood of the observed gRNA fold change for each gene differing from the distribution of the set of true neutral control gRNAs in the screen.

## Inputs: <br />
    raw counts <br />
    list of controls/ zero expressed genes (zGE) names <br />

## Module 1
Normalise raw count data as proportion per million <br />
Counts number of gRNA per gene/ intergenic <br />
Compute mean and median of normalised counts <br />

## Module 2
Compute and analyse LFC distributions of Targets and Controls (target, intergenic and zGE) <br />
Compute q-values, critical values, mean and std for zGE controls. <br />
Compute z for Targets and all CTRL, based on mean and std of zGE CTRL <br />
Determine significant gRNA enrichment/depletion accounting for both FDR and FNR rates <br />

## Running the python script
Before running the script, ensure that the config.yaml is populated.
python3 targeted_CRISPR.py   