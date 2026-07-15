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
- Created public GitHub repo `chatterjee89/4T1_hPLZF-OE_BulkRNA-seq`, pushed the
  scaffold commit.

## 2026-07-14 — Both quantification forks + downstream analysis built

- Three subagents ran in parallel and each delivered a complete unit: Fork 1
  (`forks/star_featurecounts/`, STAR+featureCounts rebuilding Plasmidsaurus's
  method), Fork 2 (`forks/salmon_kallisto/`, salmon/kallisto pseudoalignment),
  and `analysis/` (PCA/QC, per-timepoint edgeR DGE, Hallmark GSEA — run for
  real against the Plasmidsaurus matrix, the only real count source so far).
- hPLZF transgene sanity check passed cleanly: background (0–2 counts) in EV
  at every timepoint, ~0 in PLZF at D0, induced to 11–40 counts from D1
  onward (logFC 4.6–8.0, FDR < 1e-4 at D1–D3) — confirms the `sample_num`
  join, per-timepoint contrasts, and hybrid mm10+transgene reference
  reasoning are all wired correctly.
- Code review (5 finder passes across the three units) found one critical,
  confirmed bug: every `forks/star_featurecounts/` module's `conda` directive
  pointed at `${moduleDir}/envs/*.yml` instead of `${moduleDir}/../envs/*.yml`
  (envs/ is a sibling of modules/, not nested under it) — invisible under
  `-stub-run` with conda disabled, would have failed on the first real run.
  Fixed across all 12 modules. Also fixed: two analysis scripts hardcoded the
  absolute Plasmidsaurus matrix path instead of sourcing it from a gitignored
  config (added `conf/local_paths.R` + `.example`, mirroring `conf/
  local.config`'s role for R); a Cartesian-duplication edge case in
  `tximport_to_gene.R`'s gene_id/gene_name dedup; a one-directional
  join-completeness check in `build_count_matrix()`; a samplesheet
  `checkIfExists` inconsistency between the two forks.
- Also found and fixed, via direct reproduction (not just subagent reports):
  the shared `nextflow.config`'s conditional `includeConfig` pattern (both
  `def`/`if` and `try`/`catch` variants) is rejected by this Nextflow version
  (26.04.6)'s stricter config-DSL parser — switched to passing
  `conf/local.config` explicitly via `-c` on the command line. Separately,
  `${projectDir}`-based repo-root paths (`params.samplesheet`, `outdir`,
  `conda.cacheDir`) resolve to each *fork's own* directory, not the repo
  root, since each fork's main.nf is the entry script — switched to
  `${launchDir}`.
- Verified both forks end-to-end with a full, from-scratch `-stub-run`
  against the real repo config (not a workaround copy): Fork 1 196/196 stub
  tasks pass, Fork 2 76/76 pass.
- Restyled volcano plots per user request: `ggrepel::geom_text_repel` labels
  (top hits + hPLZF always), FDR+|logFC| threshold lines, square panel
  (`aspect.ratio = 1`), saved as both PNG and SVG.
- Committed in 5 logical units (scaffold; nextflow.config fix; analysis;
  Fork 1; Fork 2) and pushed all to `chatterjee89/4T1_hPLZF-OE_BulkRNA-seq`.
- Next: run Fork 1/Fork 2 for real against the full 24-sample fastq set (a
  separate, later compute job — not attempted this session), then activate
  `analysis/04_cross_method_concordance.R`'s real cross-method comparison.
