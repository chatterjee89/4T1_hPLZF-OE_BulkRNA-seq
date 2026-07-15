# analysis/R/utils.R
#
# Shared helpers for the analysis/ scripts. Designed so the *same* downstream
# code (01_sample_characterization.R, 02_dge_edgeR.R, 03_gsea_hallmark.R,
# 04_cross_method_concordance.R) can run against any of the three count-matrix
# sources this project produces:
#   - Plasmidsaurus's ground-truth STAR+featureCounts matrix (available today)
#   - forks/star_featurecounts/ output (not yet run)
#   - forks/salmon_kallisto/ output (not yet run)
#
# The only thing that differs between sources is (a) the file format/columns
# and (b) how a given sample's column(s) map back to `sample_num` in the
# sample sheet. `load_count_matrix()` below takes a `source` tag so we can
# add a branch per fork once their outputs exist, without touching any other
# script.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

#' Load a gene-level count matrix from any of this project's quantification
#' sources into a common long-form tibble:
#'   gene_id, gene_name, gene_biotype, sample_num, count
#'
#' @param path path to the source's output file
#' @param source one of "plasmidsaurus", "star_featurecounts", "salmon_kallisto"
#'
#' Only "plasmidsaurus" is implemented/testable today (Fork 1/Fork 2 haven't
#' been executed yet). The star_featurecounts / salmon_kallisto branches are
#' stubbed with a clear NotYetImplemented error and a note of the expected
#' shape, so wiring them in later is a one-function change, not a rewrite of
#' every downstream script.
load_count_matrix <- function(path, source = c("plasmidsaurus", "star_featurecounts", "salmon_kallisto")) {
  source <- match.arg(source)

  if (source == "plasmidsaurus") {
    # Wide format: gene_id, gene_name, gene_biotype, then two columns per
    # sample: BM9FTP_<sample_num>_cpm and BM9FTP_<sample_num>_count.
    # NOTE: the matrix's sample label is just BM9FTP_<sample_num> (e.g.
    # BM9FTP_10_count) -- NOT the full sample_id from the sample sheet, which
    # also carries condition/timepoint/rep (e.g. BM9FTP_10_PLZF_D1_rep1). We
    # deliberately join on sample_num later rather than sample_id string match.
    mat <- read_tsv(path, show_col_types = FALSE)

    count_cols <- grep("^BM9FTP_[0-9]+_count$", names(mat), value = TRUE)
    if (length(count_cols) == 0) {
      stop("load_count_matrix(): no *_count columns found in ", path)
    }

    long <- mat %>%
      select(gene_id, gene_name, gene_biotype, all_of(count_cols)) %>%
      pivot_longer(
        cols = all_of(count_cols),
        names_to = "sample_num",
        values_to = "count"
      ) %>%
      mutate(
        sample_num = as.integer(sub("^BM9FTP_([0-9]+)_count$", "\\1", sample_num))
      )

    return(long)
  }

  if (source == "star_featurecounts") {
    stop(
      "load_count_matrix(source = 'star_featurecounts'): not yet implemented -- ",
      "forks/star_featurecounts/ has not been executed yet, so there is no real ",
      "output to load. Expected shape once it exists: a featureCounts-style ",
      "matrix (gene_id + one raw-count column per sample, column name carrying ",
      "the sample identifier/sample_num). Implement this branch once Fork 1 ",
      "produces real output; do not fabricate placeholder data here."
    )
  }

  if (source == "salmon_kallisto") {
    stop(
      "load_count_matrix(source = 'salmon_kallisto'): not yet implemented -- ",
      "forks/salmon_kallisto/ has not been executed yet, so there is no real ",
      "output to load. Expected shape once it exists: tximport-summarised ",
      "gene-level counts (one column per sample, sample identifier/sample_num ",
      "recoverable from column names or an accompanying colData). Implement ",
      "this branch once Fork 2 produces real output; do not fabricate ",
      "placeholder data here."
    )
  }
}

#' Read the project sample sheet.
#'
#' @param path path to samplesheet/BM9FTP_samplesheet.csv
read_sample_sheet <- function(path) {
  ss <- read_csv(path, show_col_types = FALSE)
  ss$condition <- factor(ss$condition, levels = c("EV", "PLZF"))
  ss$timepoint <- factor(ss$timepoint, levels = c("D0", "D1", "D2", "D3"))
  ss
}

#' Join a long-form count tibble (from load_count_matrix()) to sample sheet
#' metadata via sample_num (NOT sample_id -- see note above), and pivot to a
#' wide gene x sample count matrix plus a matching metadata table.
#'
#' @return list(counts = matrix genes x samples, gene_info = tibble, meta = tibble)
build_count_matrix <- function(long_counts, sample_sheet) {
  joined <- long_counts %>%
    inner_join(sample_sheet, by = "sample_num")

  counts_nums <- unique(long_counts$sample_num)
  sheet_nums <- unique(sample_sheet$sample_num)
  in_counts_not_sheet <- setdiff(counts_nums, sheet_nums)
  in_sheet_not_counts <- setdiff(sheet_nums, counts_nums)
  if (length(in_counts_not_sheet) > 0) {
    warning(
      "build_count_matrix(): sample_num(s) in the count matrix with no match in ",
      "the sample sheet: ", paste(in_counts_not_sheet, collapse = ", ")
    )
  }
  if (length(in_sheet_not_counts) > 0) {
    warning(
      "build_count_matrix(): sample_num(s) in the sample sheet with no matching ",
      "count-matrix column -- these samples are silently absent from the joined ",
      "data: ", paste(in_sheet_not_counts, collapse = ", ")
    )
  }

  gene_info <- long_counts %>%
    distinct(gene_id, gene_name, gene_biotype)

  wide <- joined %>%
    select(gene_id, sample_id, count) %>%
    pivot_wider(names_from = sample_id, values_from = count)

  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$gene_id
  storage.mode(mat) <- "double"

  meta <- joined %>%
    distinct(sample_id, sample_num, condition, timepoint, replicate) %>%
    arrange(match(sample_id, colnames(mat)))

  stopifnot(identical(meta$sample_id, colnames(mat)))

  list(counts = mat, gene_info = gene_info, meta = meta)
}

#' Shared ggplot2 theme for all analysis/ plots.
theme_analysis <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90", colour = NA),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
}

#' Consistent condition color scale (EV/PLZF) reused across scripts.
condition_colors <- c(EV = "#4477AA", PLZF = "#EE6677")

#' Resolve the plots/ output directory for analysis scripts, creating it if
#' needed. Scripts are run with working directory = repo root
#' (`Rscript analysis/01_....R`), so paths are anchored there.
analysis_plots_dir <- function() {
  dir <- file.path("analysis", "plots")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir
}

#' Resolve the results/ output directory (DGE tables, GSEA tables) under analysis/.
analysis_results_dir <- function() {
  dir <- file.path("analysis", "results")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir
}

#' gene_id for the hPLZF transgene row -- used as a sanity check throughout.
HPLZF_GENE_ID <- "256578"

#' Load machine-specific data paths from conf/local_paths.R (gitignored),
#' mirroring conf/local.config's role for Nextflow but in a format R scripts
#' can source() directly. Errors with setup instructions if missing, rather
#' than falling back to a hardcoded path (which would leak an absolute local
#' path into a script committed to a public repo).
#'
#' @return an environment with the loaded variables (e.g. $PLASMIDSAURUS_MATRIX)
load_local_paths <- function() {
  path <- file.path("conf", "local_paths.R")
  if (!file.exists(path)) {
    stop(
      "conf/local_paths.R not found. Copy conf/local_paths.R.example to ",
      "conf/local_paths.R and fill in the absolute paths for this machine ",
      "(it's gitignored, so this is a one-time per-machine setup step)."
    )
  }
  env <- new.env()
  sys.source(path, envir = env)
  if (!exists("PLASMIDSAURUS_MATRIX", envir = env, inherits = FALSE)) {
    stop(
      "conf/local_paths.R exists but does not define PLASMIDSAURUS_MATRIX. ",
      "Check it against conf/local_paths.R.example -- a missing/renamed ",
      "variable here fails silently (NULL) further downstream otherwise."
    )
  }
  env
}
