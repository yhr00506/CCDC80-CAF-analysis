# Reproducibility audit — Reviewer 1 comment

Inventory only; no analysis was rerun. The archive includes existing first-party analysis scripts from project analysis, integration, validation, and revision directories, the integrated CAF metadata and recorded session information, the 19-gene list/generator, repeat summary, and per-seed status/output tables. The exact archive-to-source file mapping and SHA-256 values are in `file_manifest.json`.

Known gaps: the project root is not a Git repository and has no public remote; no DOI/public deposit currently exists. The historical R version for the three scTenifoldKnk runs is absent. Complete lockfiles for every workflow, the early CAF-state annotation source workflow, and original construction steps for the final CellChat object were not recovered in the audited paths. Some scripts contain absolute local paths. These limits preclude claiming that this candidate enables end-to-end reproduction of every manuscript result.

Checks performed: (1) source-tree inventory and SHA-256 capture; (2) CSV structural audit for integrated metadata and 19-gene file; (3) cross-check of virtual-knockout wrapper, three seed-specific status records, saved repeat summary, and generator rule.
