# 0005. Two parallel quantification forks: STAR+featureCounts and salmon/kallisto

Date: 2026-07-14
Status: Accepted

## Context

We have a ground-truth count matrix from Plasmidsaurus, produced via
STAR+featureCounts. Simply trusting it, or simply re-running the identical method,
would only tell us whether we can reproduce their pipeline — not whether the
underlying quantification is robust to methodology. We also have raw fastq, bam,
and reference material available to build our own pipelines from scratch.

## Decision

Build two independent, parallel quantification pipelines:
- `forks/star_featurecounts/`: alignment-based, mirroring Plasmidsaurus's exact
  method and tool versions (from `BM9FTP-methods.txt`) — this is a faithful
  independent reproduction, checking that we get the same answer via the same
  method.
- `forks/salmon_kallisto/`: alignment-free pseudoalignment quantification — a
  genuinely different method, so agreement here is real cross-validation, not
  just reproduction.

Both are validated against the Plasmidsaurus matrix, and both feed into the same
downstream `analysis/` scripts (PCA, edgeR, GSEA) so biological conclusions, not
just raw counts, can be compared across methods.

## Consequences

Roughly double the pipeline-building and environment-maintenance effort compared
to a single fork. In exchange, any biological finding that holds across both
pipelines is much more trustworthy than one resting on a single quantification
method.
