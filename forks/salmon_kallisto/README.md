# Fork 2: salmon / kallisto

Orthogonal, alignment-free pseudoalignment quantification, cross-validated against
both Plasmidsaurus's ground-truth matrix and Fork 1's STAR+featureCounts rebuild.
Because this is a genuinely different quantification method (not just a
reproduction of the same one), agreement with the other two sources is real
independent validation.

Pipeline: `build_transcriptome` (mm10 transcript FASTA from GENCODE vM25 + hPLZF
transcript appended, plus a `tx2gene` map) → `fastp` (trim only, **no UMI dedup** —
see `docs/adr/0006`) → `salmon index`/`salmon quant` and `kallisto index`/
`kallisto quant` in parallel → `tximport` (R) to gene-level counts.

See `docs/adr/0004`, `0005`, `0006` for the reasoning behind the reference
construction, why this fork exists alongside Fork 1, and UMI handling.

## Layout

- `main.nf` — DSL2 workflow wiring the steps below together. Reads
  `samplesheet/BM9FTP_samplesheet.csv`, ignores its `$DATA_ROOT`-placeholder
  `fastq_path` column, and instead builds real fastq paths from
  `params.fastq_dir` (from `conf/local.config`) + `sample_id` + `.fastq.gz`.
- `modules/` — one process per logical step:
  - `build_transcriptome.nf` — `gffread` extracts a spliced transcript FASTA
    from mm10 (main contigs) + GENCODE vM25, appends the hPLZF transgene as one
    extra transcript record, and emits a `tx2gene.tsv`
    (`transcript_id\tgene_id\tgene_name`) covering every mm10 transcript plus a
    synthetic `hPLZF -> hPLZF` row.
  - `fastp.nf` — trim only (polyX trim, Q20 quality trim, min length 36). No
    UMI-based dedup step anywhere in this fork (`docs/adr/0006`).
  - `salmon_index.nf` / `salmon_quant.nf` — `salmon index -k 31` (plain,
    non-decoy-aware — see note below) then `salmon quant -l A
    --validateMappings --fldMean 200 --fldSD 80` per sample.
  - `kallisto_index.nf` / `kallisto_quant.nf` — `kallisto index` then
    `kallisto quant --single --fragment-length 200 --sd 80 --plaintext` per
    sample.
  - `tximport_to_gene.nf` — runs `bin/tximport_to_gene.R`, which uses
    Bioconductor's `tximport` + the shared `tx2gene.tsv` to summarize both
    salmon's and kallisto's transcript-level estimates to gene-level counts,
    emitting two separate matrices: `salmon_gene_counts.tsv` and
    `kallisto_gene_counts.tsv` (columns: `gene_id`, `gene_name`, one column per
    sample), each joinable downstream on `gene_id` against the Plasmidsaurus
    matrix and the STAR fork's output.
- `bin/tximport_to_gene.R` — the actual tximport script (auto-added to `PATH`
  by Nextflow since it lives in this pipeline's `bin/`); reads plain relative
  filenames from its staged working directory rather than Nextflow template
  interpolation, to avoid `$`-in-R-code vs. Nextflow-template-`$` collisions.
- `envs/` — one pinned micromamba YAML per process (see versions below).

## Fragment-length assumption (salmon + kallisto)

These are single-end ~95bp reads, so neither tool can observe a true fragment
length from paired reads. Both `salmon quant` and `kallisto quant` are given
the same assumed fragment length distribution, `--fldMean/--fragment-length
200` and `--fldSD/--sd 80`, so the two tools stay comparable to each other and
not just to ground truth. **This is a documented assumption, not a measurement**
— revisit once real data is run, e.g. by cross-checking against
RSeQC/Qualimap-style metrics from the STAR fork's BAMs or vendor library-prep
documentation, and update both `salmon_quant.nf` and `kallisto_quant.nf`
together if it turns out to be off.

## Other documented simplifications

- `salmon_index.nf` builds a plain (non-decoy-aware) transcriptome index, not
  including the mm10 genome as a decoy sequence. Decoy-aware indexing is
  salmon's recommended best practice but requires indexing against the full
  ~2.7GB genome as well, a much larger one-time cost — left as a candidate
  follow-up once this fork is actually run against real data.
- Resource labels deliberately do *not* reuse the shared `alignment` label from
  `conf/base.config` (32GB/8cpu, sized for the sibling STAR fork's genome
  alignment). Pseudoalignment against a transcriptome index is much lighter,
  so these processes use `conf/base.config`'s default resources (with a modest
  `cpus 4` bump on the index/quant steps for salmon/kallisto's own
  multithreading) instead.

## Pinned tool versions

Chosen as reasonably current, stable releases as of this pipeline's build date
(2026-07-14) — these do **not** need to match Plasmidsaurus's or the STAR
fork's versions, since this fork is deliberately a different quantification
method (`docs/adr/0005`); only internal reproducibility matters here.

| Tool                    | Version | Env file             |
|-------------------------|---------|-----------------------|
| gffread                 | 0.12.7  | `envs/gffread.yml`    |
| gawk                    | 5.3.0   | `envs/gffread.yml`    |
| fastp                   | 0.23.4  | `envs/fastp.yml`      |
| salmon                  | 1.10.3  | `envs/salmon.yml`     |
| kallisto                | 0.50.1  | `envs/kallisto.yml`   |
| r-base                  | 4.3.3   | `envs/tximport.yml`   |
| bioconductor-tximport   | 1.30.0  | `envs/tximport.yml`   |
| r-readr                 | 2.1.5   | `envs/tximport.yml`   |

## Validation status

- `nextflow config forks/salmon_kallisto/main.nf` (from repo root): **fails**,
  but not because of anything in this fork — the shared top-level
  `nextflow.config`'s conditional `def localConfig = file(...); if
  (localConfig.exists()) { includeConfig ... }` pattern is incompatible with
  this machine's installed Nextflow (26.04.6)'s stricter config-DSL parser
  (`Variable declarations cannot be mixed with config statements`). Reproduced
  in complete isolation (a scratch copy with only `nextflow.config` +
  `conf/`, no fork code at all, hits the identical error), so this is a
  pre-existing, repo-wide issue in shared infrastructure outside
  `forks/salmon_kallisto/`, not something introduced here — flagging for the
  orchestrator to fix (e.g. rewrite the conditional include using the newer
  config DSL, since a bare `file(...)` call and `def` statements at the config
  top level are no longer accepted together) since it blocks `nextflow config`
  for *every* fork, not just this one.
- With that top-level config patched in an isolated scratch copy (unconditional
  `includeConfig 'conf/local.config'`, otherwise identical), `nextflow config
  forks/salmon_kallisto/main.nf` **parses cleanly** and resolves all
  `conf/local.config` params correctly.
- `nextflow run forks/salmon_kallisto/main.nf -stub-run --samplesheet
  samplesheet/BM9FTP_samplesheet.csv` (run against that same scratch copy,
  `stub:` blocks added to every process): **all 76 process instances complete
  successfully** (24 FASTP + 1 BUILD_TRANSCRIPTOME + 1 SALMON_INDEX + 1
  KALLISTO_INDEX + 24 SALMON_QUANT + 24 KALLISTO_QUANT + 1
  TXIMPORT_TO_GENE = 76), producing the expected output tree under
  `results/salmon_kallisto/{reference,fastp,salmon,kallisto,gene_counts}/`
  with correctly-named per-sample files and both `salmon_gene_counts.tsv` /
  `kallisto_gene_counts.tsv` gene-count matrices.
- Two real bugs were caught and fixed by this smoke test before it went green:
  1. Every `conda "${moduleDir}/envs/....yml"` directive was wrong —
     `moduleDir` is each module `.nf` file's own directory (`modules/`), and
     `envs/` is a *sibling* of `modules/`, not nested under it. Fixed to
     `"${moduleDir}/../envs/....yml"` throughout.
  2. `BUILD_TRANSCRIPTOME` was labeled `big_mem` (32GB) and the index/quant
     processes were labeled `alignment` (32GB/8cpu) — both borrowed from
     `conf/base.config` labels sized for the sibling STAR fork's genome-scale
     alignment. gffread/salmon/kallisto working over a transcriptome (not a
     genome) don't need that much; removed/replaced with lighter,
     right-sized resource requests (see "Other documented simplifications"
     above).
- **Important caveat about `-stub-run` itself**: on this machine's Nextflow
  (26.04.6), running `-stub-run -profile micromamba` still attempted to
  materialize a real micromamba environment for the first two processes
  reached (`gffread`, `fastp`) before hitting the memory-label error above —
  i.e. `-stub-run` does *not* skip conda/micromamba environment provisioning
  in this Nextflow version, contrary to what might be assumed. The green run
  reported above was therefore done *without* `-profile micromamba` (plain
  default profile, `conda.enabled=false`), which is sufficient to validate DAG
  wiring since every `stub:` block here is plain `bash`/`printf`/`touch` with
  no dependency on the real tools. Flagging this so the eventual real
  (non-stub) run isn't accidentally triggered by a routine stub smoke-test in
  CI or otherwise.

Nothing was run against the real 22GB+ fastq dataset or real reference files
beyond lightweight existence checks (`file(..., checkIfExists: true)` on the
real absolute paths from `conf/local.config`); no real conda/micromamba
environment was left materialized in this repo.
