#!/usr/bin/env Rscript
# analysis/03_gsea_hallmark.R
#
# Pre-ranked GSEA (fgsea) against MSigDB Hallmark gene sets, run separately
# on each timepoint's edgeR DGE result from 02_dge_edgeR.R (docs/adr/0007).
#
# Ranking metric: signed_stat = -log10(PValue) * sign(logFC).
# Rationale: this rewards genes with both a strong p-value AND a defined
# direction of change, without letting logFC magnitude alone (which can be
# noisy/inflated for low-count genes even after filterByExpr) dominate the
# ranking the way plain logFC would. Documented here per the task's ADR-style
# "your call, document it" instruction; revisit if it ever needs to match
# gseapy's default ranking exactly for cross-tool comparison.
#
# Gene sets: msigdbr's mouse-native Hallmark collection. NOTE: this msigdbr
# installation (v26.1.0) uses a newer API than the classic `category = "H"` /
# `Mm.H` convention referenced in the task brief -- mouse Hallmark is now
# fetched via db_species = "MM", collection = "MH" (see msigdbr_collections()).
# This is the direct successor of the old `Mm.H` mouse-native Hallmark set;
# functionally equivalent for this analysis's purposes.
#
# Output: one result table per timepoint under analysis/results/, and one
# dotplot of top enriched pathways per timepoint under analysis/plots/.
#
# Run from repo root:
#   Rscript analysis/03_gsea_hallmark.R [--source=plasmidsaurus]

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

source(file.path("analysis", "R", "utils.R"))

# ---- config -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
source_arg <- sub("^--source=", "", grep("^--source=", args, value = TRUE))
count_source <- if (length(source_arg)) source_arg else "plasmidsaurus"

results_dir <- analysis_results_dir()
plots_dir <- analysis_plots_dir()

cat("== 03_gsea_hallmark.R ==\n")
cat("Count source:", count_source, "\n")

# ---- Hallmark gene sets (mouse-native) -----------------------------------
hallmark_df <- msigdbr(db_species = "MM", species = "Mus musculus", collection = "MH")
hallmark_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
cat(sprintf("Loaded %d Hallmark gene sets (mouse-native, msigdbr).\n", length(hallmark_sets)))

timepoints <- c("D0", "D1", "D2", "D3")

for (tp in timepoints) {
  dge_path <- file.path(results_dir, sprintf("02_dge_%s_%s.tsv", tp, count_source))
  if (!file.exists(dge_path)) {
    cat(sprintf("\nSkipping %s: %s not found (run 02_dge_edgeR.R first).\n", tp, dge_path))
    next
  }
  cat(sprintf("\n---- Timepoint %s ----\n", tp))
  tt <- read_tsv(dge_path, show_col_types = FALSE)

  # signed ranking statistic; drop genes with no name / duplicate names by
  # keeping the most significant entry per gene symbol
  ranks_df <- tt %>%
    filter(!is.na(gene_name), gene_name != "", !is.na(PValue)) %>%
    mutate(signed_stat = -log10(pmax(PValue, .Machine$double.xmin)) * sign(logFC)) %>%
    group_by(gene_name) %>%
    slice_max(order_by = abs(signed_stat), n = 1, with_ties = FALSE) %>%
    ungroup()

  ranks <- setNames(ranks_df$signed_stat, ranks_df$gene_name)
  ranks <- sort(ranks, decreasing = TRUE)

  set.seed(42)
  gsea_res <- fgsea(pathways = hallmark_sets, stats = ranks, minSize = 10, maxSize = 500, eps = 0)
  gsea_res <- gsea_res %>% arrange(padj)

  out_path <- file.path(results_dir, sprintf("03_gsea_hallmark_%s_%s.tsv", tp, count_source))
  gsea_out <- gsea_res %>%
    mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
  write_tsv(gsea_out, out_path)
  cat("Wrote", out_path, sprintf("(%d pathways tested)\n", nrow(gsea_res)))

  n_sig <- sum(gsea_res$padj < 0.05, na.rm = TRUE)
  cat(sprintf("%d / %d Hallmark pathways significant at padj < 0.05 for %s.\n", n_sig, nrow(gsea_res), tp))

  # ---- dotplot of top enriched pathways ------------------------------------
  top_n <- 15
  plot_df <- gsea_res %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    slice_head(n = top_n) %>%
    mutate(
      pathway = sub("^HALLMARK_", "", pathway),
      pathway = reorder(pathway, NES)
    )

  if (nrow(plot_df) > 0) {
    p_dot <- ggplot(plot_df, aes(x = NES, y = pathway, size = size, color = padj)) +
      geom_point() +
      scale_color_gradient(low = "#EE6677", high = "#4477AA", name = "padj") +
      labs(
        title = sprintf("Top Hallmark GSEA pathways -- %s", tp),
        subtitle = sprintf("Source: %s | ranked by -log10(PValue)*sign(logFC) | fgsea", count_source),
        x = "Normalized Enrichment Score (NES)", y = NULL, size = "Gene set size"
      ) +
      theme_analysis()

    ggsave(file.path(plots_dir, sprintf("03_gsea_dotplot_%s_%s.png", tp, count_source)), p_dot, width = 9.5, height = 6, dpi = 150)
  } else {
    cat(sprintf("No pathways with non-NA padj for %s -- skipping dotplot.\n", tp))
  }
}

cat("\nWrote GSEA result tables to", results_dir, "and dotplots to", plots_dir, "\n")
cat("03_gsea_hallmark.R complete.\n")
