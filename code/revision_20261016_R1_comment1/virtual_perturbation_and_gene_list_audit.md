# Reviewer 1, comment 1: audit record

## 19-gene program

The project’s Supplementary Figure 3 code selects the first 19 rows after filtering the repeat summary to `gene != "CCDC80"` and `n_runs == 3`; the code does not explicitly sort before `slice_head()`. An audit found that the selected membership matches the top 19 eligible genes when independently ranked by `mean_FC`. For this revision, the export makes that ranking explicit (descending `mean_FC`, gene-symbol tie-breaker), reproducing the same membership while clarifying the rule.

Genes (ranked by mean virtual-knockout FC): **MFAP5, IGFBP6, FSTL1, FBN1, SH3BGRL3, PI16, S100A4, CLEC3B, C1QTNF3, EMP3, TPPP3, CRIP1, GSN, HTRA3, TIMP1, ELN, LINC01133, TIMP3, PTHLH.**

All 19 were present in all three computational runs, but `n_FDR05 = 0` for each in the saved repeat summary. Therefore, this is a rank/stability-selected exploratory candidate program, not a set of 19 individually FDR-significant genes. The method section and response should say this plainly.

## Virtual perturbation workflow reconstructed from the saved script and run records

- Software: `scTenifoldKnk`, version 2.1.0 (recorded in all three run-status files); the wrapper calls `scTenifoldKnk()`.
- Input: RNA counts layer from the integrated CAF Seurat object combining GSE135337, GSE190888, GSE130001, and GSE129845.
- Each run sampled 500 cells without replacement from the pooled integrated CAF object, using seeds 20260906, 20260907, and 20260908. These are computational resampling replicates from one pooled object, not independent donor-level biological replicates; cell overlap between runs is possible.
- Per run, the network input contained CCDC80 plus up to 999 Seurat variable features (1,000 genes and 500 cells in each saved run record).
- `scTenifoldKnk` arguments: `gKO = "CCDC80"`, `qc_minLSize = 0`, `nc_nNet = 5`, `nc_nCells = 100`, `nc_nComp = 3`, `td_K = 3`, `nCores = 1`, `verbose = TRUE`.
- The perturbation is an in-silico network perturbation of observational single-cell RNA-seq data. It does not establish experimental knockout effects or causal biology.

Historical runtime R and Seurat versions were not preserved in the run logs and should not be inferred from the current environment. The software version and all function parameters above are traceable to the local wrapper and per-seed status records.

## Signature comparator source audit

- CAF abundance panel: project-defined 9 genes (`COL1A1`, `COL1A2`, `DCN`, `LUM`, `COL3A1`, `COL5A1`, `SPARC`, `FAP`, `PDGFRA`).
- Narrow ECM/collagen panel: project-defined 3 genes (`COL1A1`, `COL3A1`, `FN1`).
- Broad matrisome: Matrisome Portal Homo sapiens masterlist (portal labels revision 2014; 1,027 symbols), saved under `assets/` with source manifest and checksum. This is not MatrisomeDB 2.0.
- TGFβ/myCAF: the CAF-atlas file named `TGFb_CAF` is not a specific signature: its 18-gene panel is identical across all six listed CAF states. The user selected MSigDB 2025.1 Hallmark TGFβ signaling as a pathway-level benchmark; this is not a CAF-specific myCAF signature. Do not use the repeated project panel as evidence of specificity.

## Source files

The main input hashes, Matrisome Portal source checksum, 19-gene CSV, and symbol lists are in this revision directory. Original project files have not been overwritten.
