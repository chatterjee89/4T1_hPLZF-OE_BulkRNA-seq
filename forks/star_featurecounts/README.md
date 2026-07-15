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

Status: Nextflow DSL2 pipeline implemented (`main.nf`, `modules/*.nf`,
`envs/*.yml`), `stub:` blocks added to every process, and DAG wiring smoke-
tested end-to-end with `-stub-run` (see below). **Not yet run against real
data** -- executing STAR/UMICollapse/featureCounts over the full 22GB+ fastq
set is a separate, later compute job, out of scope for this pass.

## What's implemented

`main.nf` orchestrates, per `docs/adr/0004`/`0006` and the vendor's
`methods.txt`:

1. `build_reference` (`modules/build_reference.nf`) -- concatenates
   `mm10_main.fa` + `hPLZF.fa` into one genome FASTA, measures the hPLZF
   contig length at run time (`samtools faidx`, not hardcoded) to append a
   synthetic `gene`/`transcript`/`exon` triplet (`gene_id "hPLZF"`,
   `gene_biotype "transgene"`) to a copy of the GENCODE vM25 GTF, then runs
   `STAR --runMode genomeGenerate` on the combined FASTA+GTF.
2. `gtf_to_bed12` -- converts the combined GTF to BED12 (via
   `ucsc-gtftogenepred` + `ucsc-genepredtobed`) once, for RSeQC's input.
3. `fastqc_raw` (informational) -> `fastp` (poly-X trim, 3' quality trim,
   Q15, min length 50bp, single-end).
4. `star_align` -- `--outFilterIntronMotifs RemoveNoncanonical`,
   `--outReadsUnmapped Fastx`, unsorted BAM output.
5. `samtools_sort` -- coordinate sort + index, as its own step (matching the
   vendor's method as a distinct step 4, separate from alignment).
6. `umi_dedup` -- UMICollapse `bam` mode, single-end, UMI parsed from the
   read-name suffix (`--umi-sep _`, UMICollapse's default, which matches the
   `..._TAGGCGAACACCAA` read-name layout described in the task brief).
7. `infer_strandedness` -- RSeQC `infer_experiment.py` against the dedup BAM
   + BED12; parses the `"++,--"` / `"+-,-+"` fractions into a featureCounts
   `-s` code (0/1/2). `rseqc_read_distribution` and `qualimap_rnaseq` run in
   parallel as informational-only QC (not wired into counting).
8. `featurecounts` -- **run once, batched across all 24 deduplicated BAMs**
   (see "Per-sample vs. batched featureCounts" below), `-t
   exon,three_prime_UTR -g gene_id -M --fraction`, strand flag set by a
   majority vote (in `main.nf`, no extra process) over all 24 samples'
   `infer_strandedness` results.
9. `multiqc` -- aggregates FastQC, fastp, STAR, RSeQC, Qualimap, and
   featureCounts reports.

### Per-sample vs. batched featureCounts

Chose **batched**: one `featureCounts` invocation over all 24 dedup BAMs,
producing a single `gene x sample` count matrix directly (column headers are
the staged `<sample_id>.dedup.bam` filenames). This mirrors the shape of the
Plasmidsaurus ground-truth matrix and avoids a separate matrix-merge step,
at the cost of the run only being resumable/parallelizable at the level of
"all samples" rather than per-sample. Since all 24 samples share the same
library prep, a single shared strand assignment (via majority vote) is a
reasonable simplification; if a future run shows per-sample strand
disagreement, this should be revisited (worth flagging as a follow-up if it
happens).

## Version substitutions vs. Plasmidsaurus's `methods.txt`

Checked bioconda (`https://api.anaconda.org/package/bioconda/<name>`)
directly for every tool; all versions matched exactly **except**:

- **STAR**: `methods.txt` says "v2.7.11", but bioconda only publishes the
  lettered patch builds `2.7.11a` / `2.7.11b` (no bare `2.7.11`). Substituted
  `2.7.11b` (the newer of the two, same `2.7.11` line) in both
  `envs/build_reference.yml` and `envs/star_align.yml`.
- **`ucsc-gtftogenepred` / `ucsc-genepredtobed`**: not in the vendor's tool
  list at all (this fork needs them only to feed RSeQC's BED12 requirement).
  Pinned to the latest available bioconda build (482) since there's no
  vendor version to match.

Everything else (fastqc 0.12.1, fastp 0.24.0, samtools 1.22.1, umicollapse
1.1.0, rseqc 5.0.4, qualimap 2.3, multiqc 1.32, subread 2.1.1) is available
on bioconda at the exact version in `methods.txt` and is pinned as such.

## Validation performed (this pass)

- `nextflow config forks/star_featurecounts/main.nf -profile micromamba`
  (from repo root) currently **fails**, but not because of anything in this
  fork: the shared top-level `nextflow.config` (line 9,
  `def localConfig = file(...)`) mixes a bare Groovy variable declaration
  with config statements, which Nextflow 26.04.6's config parser now rejects
  outright ("Variable declarations cannot be mixed with config statements"),
  for *every* fork, not just this one -- confirmed by running plain
  `nextflow config` from repo root with no fork argument at all, same error.
  This looks like a pre-existing incompatibility between the repo's config
  and the currently-installed Nextflow version; flagging for the
  orchestrator rather than editing `nextflow.config` myself (out of scope
  for this fork).
- Worked around it **only** for local validation: ran
  `nextflow -C <tmp-config-mirroring-nextflow.config-but-with-an-
  unconditional-includeConfig> config forks/star_featurecounts/main.nf
  -profile micromamba` against a scratch copy in `/tmp` (never touched the
  real `nextflow.config`) -- this printed a clean, fully-resolved config
  (all `params`, the `process` block, `conda.enabled/useMicromamba/cacheDir`,
  and `manifest`), confirming this fork's own Nextflow code has no config-
  level errors.
- `-stub-run`: also could not be run exactly as `-profile micromamba`,
  because Nextflow still evaluates each process's `conda` directive during
  `-stub-run` and tries to solve/build the real environment even though only
  the `stub:` block executes -- confirmed with a throwaway process/env
  pointing at a nonexistent package, which still triggered a real `libmamba`
  solve attempt. Running the *real* pinned envs under `-profile micromamba`
  would have violated the "don't materialize real conda environments" scope
  for this pass, so the smoke test was run with `conda.enabled` off (plain
  `-stub-run`, no profile) instead, using the same `/tmp` workaround config
  for the (unrelated) parsing issue above, plus a memory-label override
  (this laptop has 16GB RAM vs. `conf/base.config`'s 32GB `big_mem`/
  `alignment` labels, sized for the real run on a bigger machine).
  Result: **all 196 stub tasks completed successfully, 0 failed** --
  `build_reference` -> 24x (`fastqc_raw`, `fastp`, `star_align`,
  `samtools_sort`, `umi_dedup`, `infer_strandedness`,
  `rseqc_read_distribution`, `qualimap_rnaseq`) -> one batched
  `featurecounts` -> `multiqc`, including the strand-code majority-vote
  channel logic. Confirms the DAG wires together correctly end-to-end.
- One real bug caught and fixed by this smoke test: several `publishDir`
  directives referenced a per-sample variable (`${sample_id}`) as a plain
  interpolated string; Nextflow 26.04.6 evaluates those eagerly (unlike
  `tag`, which is always lazy) and fails with `No such variable: sample_id`.
  Fixed by wrapping all such `publishDir` values in a closure (e.g.
  `publishDir { "${params.outdir}/fastp/${sample_id}" }`) in `fastqc.nf`,
  `fastp.nf`, `star_align.nf`, `samtools_sort.nf`, `umi_dedup.nf`,
  `infer_strandedness.nf`, `rseqc_read_distribution.nf`, and `qualimap.nf`.

## Follow-up flagged for the orchestrator

- The shared `nextflow.config`'s conditional `includeConfig` pattern (line 9,
  `def localConfig = file(...)`) does not parse under the currently-installed
  Nextflow (26.04.6) and blocks `nextflow config`/`nextflow run` for *any*
  fork, not just this one. Worth a fix (and possibly an ADR note on the
  Nextflow version being pinned/required) once all forks are reviewed.
