# Original TCGA analysis addendum

This folder preserves the second-round TCGA comparison script, saved output tables, and verification records. The original script expects a local TCGA-BLCA expression object, clinical file, and project-specific paths; it is not a one-command rerun from this repository alone. R and package versions from the original run are in `sessionInfo.txt`.

For an immediately runnable check of the nine patient-level Cox comparisons using the included derived table, run `Rscript reproducibility/reproduce_tcga_models.R reproducibility` from the repository root. The derived POSTN-high myCAF comparator is a project reconstruction, not the original authors’ exact published gene list. Raw macrophage RNA-sequencing data are not included and may be requested from the corresponding author.
