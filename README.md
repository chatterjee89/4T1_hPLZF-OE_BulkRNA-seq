# 4T1_hPLZF-OE_BulkRNA-seq

Bulk RNA-seq analysis of 4T1 mouse mammary tumor cells with stable overexpression of
human ZBTB16/PLZF (hPLZF) vs. empty-vector (EV) control, sampled across 4 timepoints
(D0–D3, 3 replicates each; 24 samples total, single-end reads with UMIs).

Sequencing and an initial STAR/featureCounts analysis were performed by Plasmidsaurus;
their gene-level expression matrix is treated here as a ground-truth control. This
repo independently reproduces and cross-validates that quantification using two
parallel pipelines, then checks whether the downstream biological conclusions
(differential expression, pathway enrichment) hold up regardless of which
quantification method produced the counts.

## Repository layout

- `forks/star_featurecounts/` — Nextflow pipeline rebuilding Plasmidsaurus's method
  (fastp → STAR → UMI dedup → featureCounts) from raw fastq, pinned to their exact
  tool versions (see `docs/BM9FTP-methods.txt`-derived choices in the ADRs).
- `forks/salmon_kallisto/` — Nextflow pipeline using orthogonal pseudoalignment
  quantification (salmon + kallisto → tximport), for an independent cross-check.
- `analysis/` — sample characterization (PCA), differential expression (edgeR),
  and pathway enrichment (Hallmark GSEA), run generically against any of the three
  count-matrix sources (Plasmidsaurus / Fork 1 / Fork 2).
- `docs/adr/` — architecture decision records for the non-obvious calls made along
  the way (reference construction, UMI handling, tool choices).
- `logs/session-log.md` — running log of what was done in each work session and why.
- `conf/` — Nextflow configuration; `conf/local.config` (gitignored) holds
  machine-specific absolute paths to data/reference directories, following the
  `conf/local.config.example` template.

## What's being compared, and why

**Layer 1 — quantification concordance.** Same 24 samples, three gene-level count
matrices: Plasmidsaurus's (ground truth), our Fork 1 STAR+featureCounts rebuild
(sanity check that we reproduced their method correctly), and our Fork 2
salmon/kallisto rebuild (a genuinely different method — agreement here is real
independent validation). The `hPLZF` transgene row is checked specifically since
it's the biological readout of interest.

**Layer 2 — biological conclusion robustness.** Raw-count agreement doesn't
guarantee matching biology. Each count matrix independently goes through PCA/QC,
edgeR differential expression (EV vs PLZF per timepoint), and Hallmark GSEA — then
results are compared *across* the three sources: same genes significant, same
pathways enriched? That answers whether conclusions about what PLZF overexpression
does to these cells depend on the quantification pipeline used.

## Status

Fork 1 and Fork 2 pipelines are being built now (Nextflow code + pinned per-process
environments) but not yet run end-to-end — executing STAR/salmon over the full
fastq set is a separate, later compute job. The `analysis/` scripts are written
generically and proven out today against the one real count matrix that already
exists (Plasmidsaurus's).

## Reference & data

- Genome: UCSC mm10 (`mm10_main.fa`, 21 main contigs) + GENCODE vM25 annotation,
  matching what Plasmidsaurus used.
- A custom `hPLZF` transgene contig (human ZBTB16 cDNA) is appended to the
  reference so reads from the overexpression construct are captured rather than
  lost or misassigned to the endogenous mouse `Zbtb16` locus.
- Raw fastq/bam and reference files are large (tens of GB) and live outside this
  repo; see `conf/local.config.example` for how paths are wired in.
