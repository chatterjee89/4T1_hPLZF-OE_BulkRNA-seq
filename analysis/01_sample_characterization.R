#!/usr/bin/env Rscript
# analysis/01_sample_characterization.R
#
# Sample-level QC gate: TMM-normalized PCA + sample-sample correlation heatmap,
# colored by condition (EV/PLZF) and faceted/shaped by timepoint (D0-D3).
#
# Run from repo root:
#   Rscript analysis/01_sample_characterization.R [--source=plasmidsaurus]
#
# Generic across count-matrix sources: which source to load is the only thing
# that changes (see analysis/R/utils.R::load_count_matrix()). Today only the
# Plasmidsaurus matrix exists; Fork 1/Fork 2 aren't run yet.

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
})

source(file.path("analysis", "R", "utils.R"))

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

cat("== 01_sample_characterization.R ==\n")
cat("Count source:", count_source, "\n")

# ---- load data ----------------------------------------------------------
long_counts <- load_count_matrix(matrix_paths[[count_source]], source = count_source)
sample_sheet <- read_sample_sheet(sample_sheet_path)
built <- build_count_matrix(long_counts, sample_sheet)

counts <- built$counts
gene_info <- built$gene_info
meta <- built$meta

cat(sprintf("Loaded %d genes x %d samples.\n", nrow(counts), ncol(counts)))

# ---- DGEList + TMM normalization ----------------------------------------
dge <- DGEList(counts = counts, samples = meta)
keep <- filterByExpr(dge, group = meta$condition)
cat(sprintf("filterByExpr: keeping %d / %d genes.\n", sum(keep), length(keep)))
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# ---- PCA on log-CPM -------------------------------------------------------
pca <- prcomp(t(logcpm), center = TRUE, scale. = FALSE)
var_explained <- (pca$sdev^2) / sum(pca$sdev^2) * 100

pca_df <- as.data.frame(pca$x[, 1:min(4, ncol(pca$x))]) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(meta, by = "sample_id")

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition, shape = timepoint)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = condition_colors) +
  labs(
    title = "PCA of log-CPM expression (TMM-normalized)",
    subtitle = sprintf("Source: %s", count_source),
    x = sprintf("PC1 (%.1f%% var)", var_explained[1]),
    y = sprintf("PC2 (%.1f%% var)", var_explained[2]),
    color = "Condition", shape = "Timepoint"
  ) +
  theme_analysis()

ggsave(file.path(plots_dir, sprintf("01_pca_%s.png", count_source)), p_pca, width = 7, height = 5.5, dpi = 150)

p_pca_facet <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = condition_colors) +
  facet_wrap(~timepoint) +
  labs(
    title = "PCA of log-CPM expression, faceted by timepoint",
    subtitle = sprintf("Source: %s", count_source),
    x = sprintf("PC1 (%.1f%% var)", var_explained[1]),
    y = sprintf("PC2 (%.1f%% var)", var_explained[2]),
    color = "Condition"
  ) +
  theme_analysis()

ggsave(file.path(plots_dir, sprintf("01_pca_facet_timepoint_%s.png", count_source)), p_pca_facet, width = 8, height = 6.5, dpi = 150)

# ---- sample-sample correlation heatmap -----------------------------------
sample_cor <- cor(logcpm, method = "pearson")

annot <- meta %>%
  select(sample_id, condition, timepoint) %>%
  tibble::column_to_rownames("sample_id")

annot_colors <- list(
  condition = condition_colors,
  timepoint = setNames(
    scales::hue_pal()(nlevels(meta$timepoint)),
    levels(meta$timepoint)
  )
)

png(file.path(plots_dir, sprintf("01_sample_correlation_heatmap_%s.png", count_source)),
    width = 1400, height = 1300, res = 150)
pheatmap(
  sample_cor,
  annotation_col = annot,
  annotation_colors = annot_colors,
  main = sprintf("Sample-sample Pearson correlation (log-CPM), source: %s", count_source),
  fontsize = 8
)
dev.off()

# ---- hPLZF transgene sanity check (raw view, quick look) ------------------
if (HPLZF_GENE_ID %in% rownames(dge$counts)) {
  hplzf_cpm <- cpm(dge)[HPLZF_GENE_ID, ]
  hplzf_df <- data.frame(sample_id = names(hplzf_cpm), cpm = hplzf_cpm) %>%
    left_join(meta, by = "sample_id")
  cat("\n--- hPLZF (gene_id 256578) CPM by condition (quick look; see 02_dge_edgeR.R for stats) ---\n")
  print(hplzf_df %>% group_by(condition) %>% summarise(mean_cpm = mean(cpm), .groups = "drop"))
} else {
  cat(
    "\nNote: hPLZF (gene_id 256578) was dropped by the genome-wide filterByExpr(group = condition)\n",
    "call used for this PCA/QC pass. This is expected, not a bug: with 12 EV vs 12 PLZF samples\n",
    "(> filterByExpr's large.n = 10), the default min.prop = 0.7 inflates the required number of\n",
    "expressing samples to ceil(10 + (12-10)*0.7) = 12/12 -- i.e. it demands expression in *every*\n",
    "sample of a group. hPLZF is only induced from D1 onward (~0 at D0 even in PLZF samples; see\n",
    "raw counts below), so it fails this genome-wide filter. Per-timepoint contrasts in\n",
    "02_dge_edgeR.R use 3-vs-3 sample groups (below large.n), where this inflation does not apply,\n",
    "so hPLZF is expected to pass filtering there for D1/D2/D3 (and correctly drop out at D0).\n"
  )
}

cat("\nWrote plots to", plots_dir, "\n")
cat("01_sample_characterization.R complete.\n")
