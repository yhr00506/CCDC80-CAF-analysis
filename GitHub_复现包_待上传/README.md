# CCDC80-associated CAF analysis: code and reproducibility materials

This repository accompanies *Exploratory CCDC80-Associated CAF Extracellular-Matrix Program and Myeloid-Context Associations in Bladder Cancer* (Biomedicines manuscript biomedicines-4606757). The upload target is https://github.com/yhr00506/CCDC80-CAF-analysis. It contains archived analysis scripts, integrated CAF metadata, the exploratory 19-gene selection records, and a portable TCGA patient-level model rerun.

## Contents

- `code/`: archived first-party analysis scripts. Some retain project-specific absolute paths and require configuration before reuse.
- `metadata/`: integrated CAF candidate metadata and available session information.
- `signature/`: three saved CCDC80 virtual-knockout resampling runs, repeat summary, selection script, and final 19-gene list.
- `reproducibility/`: portable 19-gene selection check, 400-patient TCGA model rerun, input table, saved comparison values, and outputs.
- `provenance/`: source-file inventory and known historical workflow gaps.

## Run the portable checks

From the repository root:

```bash
Rscript reproducibility/reproduce_signature_selection.R .
Rscript reproducibility/reproduce_tcga_models.R reproducibility
```

The first command verifies the 19 selected symbols against the saved repeat summary. The second independently fits nine nested TCGA Cox comparisons from the supplied derived patient-level table and checks their likelihood-ratio statistics and adjusted P values against the saved results. The latter requires the R `survival` package; the original analysis used R 4.5.1 and `survival` 3.8-3.

## Scope and data availability

The 19 genes are exploratory and rank/stability-selected. None met FDR < 0.05 in the saved virtual-knockout runs; the three runs are computational resamples, not donor-level replicates.

The portable TCGA rerun begins with a derived 400-patient table. Its upstream construction used the TCGA-BLCA `unstranded` count assay, log₂(count + 1) transformation without CPM/TPM or library-size normalization, within-patient averaging of duplicate primary-tumor aliquots, gene-wise z standardization across 406 unique patient-level profiles, and score standardization in the 400-patient survival set. The original expression extraction script and saved outputs are in `reproducibility/tcga_addendum/`; that script needs source data and local path configuration.

This repository does not include the raw TCGA expression object, raw experimental macrophage RNA-seq data, or a complete one-command rerun of every historical CAF-state and CellChat analysis. Public source accessions are listed in the manuscript. Experimental macrophage RNA-seq data may be requested from the corresponding author. See `provenance/reproducibility_audit.md` for the historical workflow boundary.
