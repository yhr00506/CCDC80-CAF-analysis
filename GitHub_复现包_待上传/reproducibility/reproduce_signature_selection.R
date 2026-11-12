arguments <- commandArgs(trailingOnly = TRUE)
repository_root <- if (length(arguments) >= 1) arguments[[1]] else "."
summary_file <- file.path(repository_root, "signature/analysis_revisions_20260906/virtual_knockout_repeat/CCDC80/CCDC80_repeat_summary.csv")
published_file <- file.path(repository_root, "signature/revision_20261016_R1_comment1/assets/CCDC80_associated_19_gene_symbols.txt")
stopifnot(file.exists(summary_file), file.exists(published_file))

summary <- read.csv(summary_file, check.names = FALSE)
candidates <- summary[summary$gene != "CCDC80" & summary$n_runs == 3, , drop = FALSE]
candidates <- candidates[order(-candidates$mean_FC, candidates$gene), , drop = FALSE]
selected <- head(candidates, 19)
expected <- readLines(published_file)
stopifnot(nrow(selected) == 19, identical(as.character(selected$gene), expected))
stopifnot(all(selected$n_FDR05 == 0))
selected$rank <- seq_len(nrow(selected))
write.csv(selected, file.path(repository_root, "reproducibility/selected_19_genes_verified.csv"), row.names = FALSE)
cat("Verified all 19 rank/stability-selected genes; none has FDR < 0.05 in the saved runs.\n")
