# Session log

## 2026-07-14 — Project bootstrap

- Explored existing data on disk (`~/Documents/Work/Computational/BulkRNAseq/`):
  24 fastq (single-end, UMI in read name), Plasmidsaurus BAMs + expression matrix,
  mm10 + GENCODE vM25 reference, hPLZF transgene FASTA, methods.txt describing
  Plasmidsaurus's exact pipeline (fastp → STAR → UMIcollapse → featureCounts →
  edgeR → gseapy Hallmark GSEA).
- Decided project architecture: Nextflow orchestration, micromamba per-process
  envs, two parallel quantification forks (STAR+featureCounts rebuilding
  Plasmidsaurus's method; salmon/kallisto as an orthogonal method), plus a shared
  `analysis/` stage (PCA, edgeR DGE, Hallmark GSEA) applied identically to all
  three count-matrix sources. See `docs/adr/0001`–`0007`.
- Installed toolchain on this Mac: Homebrew, OpenJDK (via `JAVA_HOME`, no sudo
  symlink needed), Nextflow 26.04.6, `gh` CLI (authenticated as `chatterjee89`),
  micromamba.
- Scaffolded repo structure (this commit): README, CLAUDE.md, ADRs, this log,
  Nextflow configs, sample sheet, empty fork/analysis directories with READMEs,
  `.claude/skills/new-analytical-fork`.
- Next: three subagents build `forks/star_featurecounts/`, `forks/salmon_kallisto/`,
  and `analysis/` (the latter run now against the real Plasmidsaurus matrix, since
  the quantification forks aren't executed yet — that's a separate later compute
  job).
