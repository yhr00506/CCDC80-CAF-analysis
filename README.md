# Supplementary Code and Metadata 1 — local release candidate

This archive is a local review copy prepared from project files. It has not been uploaded to GitHub, Zenodo, or another public repository. The authors plan to deposit the finalized archive in Zenodo before publication and add its DOI/URL to the manuscript.

## Contents
- `code/`: first-party analysis scripts found in the project analysis/revision directories; original source files are retained. Some scripts refer to project-specific absolute paths and require path configuration before reuse.
- `metadata/`: integrated CAF-specific cell metadata (5,589 rows, 31 columns), candidate audit, and recorded session information.
- `signature/`: 19-gene program table and symbol list, generation script/audit, repeat summary and per-seed status/differential-regulation outputs.
- `provenance/`: SHA-256 inventory and scope/known-gap statement.

The 19 symbols are MFAP5, IGFBP6, FSTL1, FBN1, SH3BGRL3, PI16, S100A4, CLEC3B, C1QTNF3, EMP3, TPPP3, CRIP1, GSN, HTRA3, TIMP1, ELN, LINC01133, TIMP3, and PTHLH. The selection uses genes present in all three computational resampling runs, ranked by mean virtual-knockout fold change (descending; symbol ascending for ties), excluding CCDC80; none of the 19 had FDR < 0.05 in the saved repeat summary. The runs use 500 sampled cells each and are not donor-level biological replicates.

## Known reproducibility boundaries
The audit did not recover a complete historical environment lockfile, the R runtime version for the virtual-knockout runs, all early CAF-state annotation steps, or original construction steps for the final CellChat object. Full end-to-end reproduction therefore remains pending recovery of these records. Required input expression objects are not redistributed; the source public accessions are documented in the manuscript. Raw expression data and large RDS perturbation objects are omitted. Please consult `provenance/reproducibility_audit.md` before interpreting this package as a complete rerun workflow.

No license is asserted for this candidate archive; licensing must be selected by the authors before public deposit.
