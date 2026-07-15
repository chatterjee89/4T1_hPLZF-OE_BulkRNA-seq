#!/usr/bin/env Rscript
#
# Summarizes salmon's and kallisto's transcript-level quantifications to
# gene-level counts using the shared tx2gene map (mm10 GENCODE vM25 + hPLZF),
# via Bioconductor's tximport. Produces two separate gene-count matrices
# (one per quantifier) so each can be joined later against the Plasmidsaurus
# ground-truth matrix and the STAR/featureCounts fork's output on gene_id.
#
# Invoked with no arguments by modules/tximport_to_gene.nf. Nextflow's `conda`
# process directive stages this pipeline's bin/ on PATH, and the process's
# `path` inputs are staged directly into the task's working directory (which
# is this script's cwd when it runs), so everything below is read via plain
# relative paths:
#
#   tx2gene.tsv           <- BUILD_TRANSCRIPTOME.out.tx2gene
#   <sample_id>_salmon/quant.sf       <- one per sample, from SALMON_QUANT
#   <sample_id>_kallisto/abundance.tsv <- one per sample, from KALLISTO_QUANT

suppressMessages({
  library(tximport)
  library(readr)
})

tx2gene_full <- read_tsv(
  "tx2gene.tsv",
  col_names = c("transcript_id", "gene_id", "gene_name"),
  col_types = "ccc"
)
tx2gene    <- as.data.frame(tx2gene_full[, c("transcript_id", "gene_id")])
# Keep the FIRST gene_name seen per gene_id, not unique(gene_id, gene_name)
# pairs -- the latter would silently Cartesian-duplicate a gene's output row
# (same count, listed twice) if a gene_id ever maps to >1 gene_name spelling
# in the GTF (not the case in the current mm10 vM25 GTF, but not guaranteed
# to stay that way).
gene_names <- as.data.frame(tx2gene_full[, c("gene_id", "gene_name")])
gene_names <- gene_names[!duplicated(gene_names$gene_id), ]

salmon_dirs <- sort(list.dirs(".", recursive = FALSE, full.names = FALSE))
salmon_dirs <- salmon_dirs[grepl("_salmon$", salmon_dirs)]
if (length(salmon_dirs) == 0) stop("No *_salmon quant directories found in cwd.")
salmon_samples <- sub("_salmon$", "", salmon_dirs)
salmon_files <- file.path(salmon_dirs, "quant.sf")
names(salmon_files) <- salmon_samples

kallisto_dirs <- sort(list.dirs(".", recursive = FALSE, full.names = FALSE))
kallisto_dirs <- kallisto_dirs[grepl("_kallisto$", kallisto_dirs)]
if (length(kallisto_dirs) == 0) stop("No *_kallisto quant directories found in cwd.")
kallisto_samples <- sub("_kallisto$", "", kallisto_dirs)
kallisto_files <- file.path(kallisto_dirs, "abundance.tsv")
names(kallisto_files) <- kallisto_samples

# Writes a tximport `counts` matrix (transcript-level estimates already
# summed to gene level by tximport itself, via tx2gene) out as a TSV with
# gene_id + gene_name leading columns, one sample column per fastq, sorted
# by gene_id -- matching the shape of the Plasmidsaurus ground-truth matrix
# (gene_id, gene_name, gene_biotype, <samples...>) closely enough to join on
# gene_id downstream.
write_gene_matrix <- function(txi, out_path) {
  counts <- as.data.frame(round(txi$counts))
  counts$gene_id <- rownames(counts)
  counts <- merge(gene_names, counts, by = "gene_id", all.y = TRUE)
  sample_cols <- setdiff(colnames(counts), c("gene_id", "gene_name"))
  counts <- counts[order(counts$gene_id), c("gene_id", "gene_name", sample_cols)]
  write.table(counts, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

txi_salmon <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write_gene_matrix(txi_salmon, "salmon_gene_counts.tsv")

txi_kallisto <- tximport(kallisto_files, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
write_gene_matrix(txi_kallisto, "kallisto_gene_counts.tsv")
