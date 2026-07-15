# Fork 1: STAR + featureCounts

Rebuilds Plasmidsaurus's quantification method from raw fastq (not their BAMs),
pinned to their exact tool versions, so the result is an independently reproduced
count matrix comparable to their ground-truth expression matrix.

Pipeline: `build_reference` (mm10 + GENCODE vM25 + hPLZF transgene contig) →
`fastp` v0.24.0 → `STAR` v2.7.11 → `samtools` v1.22.1 sort → `UMIcollapse` v1.1.0
dedup → `infer_strandedness` (RSeQC) → `featureCounts` (subread v2.1.1) →
`MultiQC` v1.32.

See `docs/adr/0004`, `0005`, `0006` for the reasoning behind the reference
construction, why this fork exists alongside Fork 2, and UMI handling.

Status: scaffolding in progress (`main.nf` / `modules/` / `envs/` to follow).
