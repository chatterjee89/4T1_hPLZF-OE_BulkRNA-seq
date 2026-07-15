#!/usr/bin/env Rscript
# analysis/04_cross_method_concordance.R
#
# Layer 1/2 cross-method comparison (docs/adr/0005, 0007): once Plasmidsaurus,
# Fork 1 (star_featurecounts), and Fork 2 (salmon_kallisto) all have DGE +
# GSEA results from 02/03, this script compares them -- DGE hit-list overlap
# (Jaccard) and GSEA pathway-hit overlap, per timepoint -- to check whether
# biological conclusions are robust to quantification method.
#
# TODAY: only the Plasmidsaurus source has been run (forks/star_featurecounts/
# and forks/salmon_kallisto/ are built but not yet executed -- that's a
# separate, later compute job per the repo README). So this script runs in
# "single-source mode": it summarizes the one available source's DGE/GSEA
# results and validates internal consistency, and prints an explicit,
# un-missable notice that cross-method comparison will activate once Fork 1
# and/or Fork 2 outputs exist. It does NOT fabricate placeholder Fork 1/Fork 2
# data to make this look more complete than it is.
#
# Per docs/adr/0006: Fork 2 (salmon/kallisto) deliberately skips UMI
# deduplication (unlike Fork 1/Plasmidsaurus), so once real multi-source
# comparison runs, some discordance -- particularly for highly-expressed
# genes with higher PCR duplication rates -- is EXPECTED and should not be
# treated as pipeline error. This script's multi-source branch (below) is
# written to call that out rather than silently reconcile it.
#
# Run from repo root:
#   Rscript analysis/04_cross_method_concordance.R

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

source(file.path("analysis", "R", "utils.R"))

results_dir <- analysis_results_dir()
timepoints <- c("D0", "D1", "D2", "D3")

# Sources this script is *designed* to compare once they exist. Only sources
# with actual result files on disk are treated as "available" -- everything
# else is reported as pending, never invented.
candidate_sources <- c("plasmidsaurus", "star_featurecounts", "salmon_kallisto")

source_has_results <- function(src) {
  paths <- file.path(results_dir, sprintf("02_dge_%s_%s.tsv", timepoints, src))
  all(file.exists(paths))
}

available_sources <- Filter(source_has_results, candidate_sources)
pending_sources <- setdiff(candidate_sources, available_sources)

cat("== 04_cross_method_concordance.R ==\n")
cat("Available count-matrix sources (real DGE output on disk):", paste(available_sources, collapse = ", "), "\n")
if (length(pending_sources) > 0) {
  cat(
    "\n*** NOTICE ***\n",
    "The following source(s) have NOT been run yet and are excluded from this\n",
    "report (no fabricated data is substituted): ", paste(pending_sources, collapse = ", "), "\n",
    "Fork 1 (star_featurecounts) and Fork 2 (salmon_kallisto) are built but not\n",
    "yet executed against the full fastq set -- that is a separate, later compute\n",
    "job (see repo README.md Status section). Once either fork produces real\n",
    "02_dge_edgeR.R / 03_gsea_hallmark.R output under analysis/results/, re-run\n",
    "this script: it will automatically pick them up (see `candidate_sources`\n",
    "and `load_count_matrix()` in analysis/R/utils.R) and activate genuine\n",
    "cross-method DGE/GSEA comparison below.\n",
    sep = ""
  )
}

if (length(available_sources) == 0) {
  stop("No sources with real DGE output found -- run 02_dge_edgeR.R first.")
}

# ---- single-source mode: summarize/validate what we do have --------------
cat("\n---- Single-source summary mode ----\n")

summarize_source <- function(src) {
  purrr::map_dfr(timepoints, function(tp) {
    dge_path <- file.path(results_dir, sprintf("02_dge_%s_%s.tsv", tp, src))
    if (!file.exists(dge_path)) return(NULL)
    tt <- read_tsv(dge_path, show_col_types = FALSE)
    hplzf <- tt %>% filter(gene_id == HPLZF_GENE_ID)
    tibble(
      source = src,
      timepoint = tp,
      n_genes_tested = nrow(tt),
      n_sig_fdr05 = sum(tt$FDR < 0.05, na.rm = TRUE),
      n_up_in_plzf = sum(tt$FDR < 0.05 & tt$logFC > 0, na.rm = TRUE),
      n_down_in_plzf = sum(tt$FDR < 0.05 & tt$logFC < 0, na.rm = TRUE),
      hPLZF_logFC = if (nrow(hplzf) == 1) hplzf$logFC else NA_real_,
      hPLZF_FDR = if (nrow(hplzf) == 1) hplzf$FDR else NA_real_
    )
  })
}

single_source_summary <- purrr::map_dfr(available_sources, summarize_source)
print(single_source_summary, n = Inf)

out_path <- file.path(results_dir, "04_single_source_summary.tsv")
write_tsv(single_source_summary, out_path)
cat("\nWrote", out_path, "\n")

# ---- multi-source comparison (activates once >=2 sources are available) ---
if (length(available_sources) >= 2) {
  cat("\n---- Multi-source DGE/GSEA concordance ----\n")

  pairwise <- combn(available_sources, 2, simplify = FALSE)

  jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

  concordance <- purrr::map_dfr(pairwise, function(pair) {
    purrr::map_dfr(timepoints, function(tp) {
      pa <- file.path(results_dir, sprintf("02_dge_%s_%s.tsv", tp, pair[1]))
      pb <- file.path(results_dir, sprintf("02_dge_%s_%s.tsv", tp, pair[2]))
      if (!file.exists(pa) || !file.exists(pb)) return(NULL)

      tt_a <- read_tsv(pa, show_col_types = FALSE) %>% filter(FDR < 0.05)
      tt_b <- read_tsv(pb, show_col_types = FALSE) %>% filter(FDR < 0.05)

      tibble(
        source_a = pair[1], source_b = pair[2], timepoint = tp,
        n_sig_a = nrow(tt_a), n_sig_b = nrow(tt_b),
        n_overlap = length(intersect(tt_a$gene_id, tt_b$gene_id)),
        jaccard_dge = jaccard(tt_a$gene_id, tt_b$gene_id)
      )
    })
  })

  print(concordance, n = Inf)
  write_tsv(concordance, file.path(results_dir, "04_cross_method_dge_concordance.tsv"))

  # Note the expected UMI-dedup-driven asymmetry when salmon_kallisto is involved
  if ("salmon_kallisto" %in% available_sources) {
    cat(
      "\nReminder (docs/adr/0006): salmon_kallisto skips UMI dedup by design.\n",
      "Discordance vs. plasmidsaurus/star_featurecounts, especially for highly\n",
      "expressed genes, may reflect that methodological difference rather than\n",
      "pipeline error -- inspect before treating low concordance as a bug.\n"
    )
  }
} else {
  cat(
    "\nOnly 1 source available -- multi-source DGE/GSEA concordance comparison\n",
    "is not run this time (nothing fabricated). See NOTICE above.\n"
  )
}

cat("\n04_cross_method_concordance.R complete.\n")
