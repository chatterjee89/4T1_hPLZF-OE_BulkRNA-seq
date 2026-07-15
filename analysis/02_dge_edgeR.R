#!/usr/bin/env Rscript
# analysis/02_dge_edgeR.R
#
# edgeR differential expression, EV vs PLZF, fit SEPARATELY per timepoint
# (D0, D1, D2, D3 -- 4 independent 3-vs-3 contrasts, not one combined model).
# TMM normalization + filterByExpr defaults, matching Plasmidsaurus's own
# method (docs/adr/0007) for direct comparability.
#
# Test used: glmQLFit / glmQLFTest (quasi-likelihood F-test), not exactTest.
# Rationale: glmQLFit's empirical Bayes moderation of the QL dispersion is
# more robust with small (n=3 per group) sample sizes than exactTest, which
# assumes the tagwise dispersion estimate is exact; QL also gives a
# straightforward path to more complex designs later (e.g. adding covariates)
# without switching frameworks. Recommended as the default modern edgeR
# workflow.
#
# Output: one result table per timepoint (gene_id, gene_name, gene_biotype,
# logFC, logCPM, F, PValue, FDR) under analysis/results/, and one volcano
# plot per timepoint under analysis/plots/. The hPLZF transgene row
# (gene_id 256578) is printed explicitly per timepoint as a sanity check.
#
# Run from repo root:
#   Rscript analysis/02_dge_edgeR.R [--source=plasmidsaurus]

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(svglite)
  library(dplyr)
  library(readr)
})

source(file.path("analysis", "R", "utils.R"))

# Volcano plot thresholds/style (see "volcano plot" section below).
LFC_THRESH <- 1
FDR_THRESH <- 0.05

# ---- config -----------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
source_arg <- sub("^--source=", "", grep("^--source=", args, value = TRUE))
count_source <- if (length(source_arg)) source_arg else "plasmidsaurus"

local_paths <- load_local_paths()
matrix_paths <- list(
  plasmidsaurus = local_paths$PLASMIDSAURUS_MATRIX
)
sample_sheet_path <- file.path("samplesheet", "BM9FTP_samplesheet.csv")

plots_dir <- analysis_plots_dir()
results_dir <- analysis_results_dir()

cat("== 02_dge_edgeR.R ==\n")
cat("Count source:", count_source, "\n")

# ---- load data ----------------------------------------------------------
long_counts <- load_count_matrix(matrix_paths[[count_source]], source = count_source)
sample_sheet <- read_sample_sheet(sample_sheet_path)
built <- build_count_matrix(long_counts, sample_sheet)

counts <- built$counts
gene_info <- built$gene_info
meta <- built$meta

timepoints <- levels(meta$timepoint)

hplzf_summary <- list()
all_results <- list()

for (tp in timepoints) {
  cat(sprintf("\n---- Timepoint %s ----\n", tp))
  idx <- meta$timepoint == tp
  sub_counts <- counts[, idx, drop = FALSE]
  sub_meta <- meta[idx, , drop = FALSE]
  sub_meta$condition <- droplevels(sub_meta$condition)

  cond_tab <- table(sub_meta$condition)
  cat(sprintf("n = %d samples (%d EV / %d PLZF)\n", ncol(sub_counts), cond_tab["EV"], cond_tab["PLZF"]))

  dge <- DGEList(counts = sub_counts, samples = sub_meta, genes = gene_info[match(rownames(sub_counts), gene_info$gene_id), ])
  keep <- filterByExpr(dge, group = sub_meta$condition)
  cat(sprintf("filterByExpr: keeping %d / %d genes for %s.\n", sum(keep), length(keep), tp))
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")

  design <- model.matrix(~condition, data = sub_meta)
  colnames(design) <- make.names(colnames(design))

  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = ncol(design))  # PLZF vs EV (EV is reference level)

  tt <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  tt <- tt %>%
    select(gene_id, gene_name, gene_biotype, logFC, logCPM, F, PValue, FDR) %>%
    arrange(PValue)

  out_path <- file.path(results_dir, sprintf("02_dge_%s_%s.tsv", tp, count_source))
  write_tsv(tt, out_path)
  cat("Wrote", out_path, sprintf("(%d genes tested)\n", nrow(tt)))

  all_results[[tp]] <- tt

  # ---- hPLZF sanity check -------------------------------------------------
  hplzf_row <- tt %>% filter(gene_id == HPLZF_GENE_ID)
  if (nrow(hplzf_row) >= 1) {
    if (nrow(hplzf_row) > 1) {
      warning(sprintf(
        "hPLZF sanity check [%s]: gene_id %s appears %d times in the result table (expected 1) -- reporting the first row; this points at a data-quality issue upstream, not a normal 'filtered out' case.",
        tp, HPLZF_GENE_ID, nrow(hplzf_row)
      ))
      hplzf_row <- hplzf_row[1, ]
    }
    cat(sprintf(
      "hPLZF sanity check [%s]: logFC = %.3f, PValue = %.3g, FDR = %.3g\n",
      tp, hplzf_row$logFC, hplzf_row$PValue, hplzf_row$FDR
    ))
    hplzf_summary[[tp]] <- data.frame(timepoint = tp, logFC = hplzf_row$logFC, PValue = hplzf_row$PValue, FDR = hplzf_row$FDR)
  } else {
    cat(sprintf(
      "hPLZF sanity check [%s]: gene_id %s was filtered out by filterByExpr at this timepoint\n",
      tp, HPLZF_GENE_ID
    ))
    hplzf_summary[[tp]] <- data.frame(timepoint = tp, logFC = NA, PValue = NA, FDR = NA)
  }

  # ---- volcano plot --------------------------------------------------------
  # House style: color by Sig (FDR + |logFC| threshold), label top hits plus
  # hPLZF always (regardless of significance, since it's the built-in sanity
  # check), square panel, saved as both PNG and SVG.
  d <- tt %>%
    mutate(Sig = case_when(
      FDR < FDR_THRESH & logFC > LFC_THRESH ~ "Up",
      FDR < FDR_THRESH & logFC < -LFC_THRESH ~ "Down",
      TRUE ~ "NS"
    ))

  n_up <- sum(d$Sig == "Up")
  n_down <- sum(d$Sig == "Down")
  y_max <- max(-log10(d$PValue), na.rm = TRUE) * 1.05

  pal <- c(Up = "#EE6677", Down = "#4477AA", NS = "grey70")

  lbl_df <- bind_rows(
    d %>% filter(Sig != "NS") %>% arrange(PValue) %>% head(30),
    d %>% filter(gene_id == HPLZF_GENE_ID)
  ) %>% distinct(gene_id, .keep_all = TRUE)

  base_size <- 11
  day <- tp

  p_volcano <- ggplot(d, aes(x = logFC, y = -log10(PValue), color = Sig)) +
    geom_point(size = 0.55, alpha = 0.65) +
    geom_text_repel(
      data = lbl_df, aes(label = gene_name),
      size = 1.9, segment.size = 0.15, max.overlaps = 60,
      box.padding = 0.2, point.padding = 0.05,
      segment.alpha = 0.4, show.legend = FALSE) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH),
               linetype = "dashed", color = "gray50", linewidth = 0.3) +
    geom_hline(yintercept = -log10(FDR_THRESH),
               linetype = "dashed", color = "gray50", linewidth = 0.3) +
    scale_color_manual(values = pal) +
    scale_y_continuous(limits = c(0, y_max)) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 2.3,
             label = sprintf("Up: %d   Down: %d", n_up, n_down)) +
    labs(
      title = sprintf("4T1 %s -- PLZF vs EV", day),
      x = "log2 Fold Change",
      y = "-log10(p-value)") +
    theme_bw(base_size = base_size) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = base_size),
      legend.position = "none",
      panel.grid      = element_blank(),
      aspect.ratio    = 1,  # square inside plot area, independent of data x/y ranges
      plot.margin     = margin(4, 4, 4, 4))

  for (ext in c("png", "svg")) {
    ggsave(
      file.path(plots_dir, sprintf("02_volcano_%s_%s.%s", tp, count_source, ext)),
      p_volcano, width = 5, height = 5, dpi = 150
    )
  }
}

cat("\n==== hPLZF (gene_id 256578) sanity-check summary across timepoints ====\n")
hplzf_all <- bind_rows(hplzf_summary)
print(hplzf_all, row.names = FALSE)
write_tsv(hplzf_all, file.path(results_dir, sprintf("02_hplzf_sanity_check_%s.tsv", count_source)))

cat("\nWrote result tables to", results_dir, "and volcano plots to", plots_dir, "\n")
cat("02_dge_edgeR.R complete.\n")
