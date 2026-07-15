# Working conventions for this repo

## Workflow

For any non-trivial change: **explore → plan → implement → commit**.
- Explore the relevant code/data before proposing anything.
- Always pause and ask clarifying questions (interview) before finalizing a plan —
  do not assume requirements, especially anything touching reference construction,
  UMI/dedup handling, or statistical method choices.
- Implement, then show the diff before committing.
- Commit in small, reviewable units — one fork/module/analysis stage per commit,
  not one giant commit for the whole repo.

## Data & reference handling

- Never commit raw data (fastq/bam), reference genomes/annotations, or anything
  matched by `.gitignore`. These live outside the repo and are referenced via
  `conf/local.config` (gitignored; see `conf/local.config.example` for the
  template that new machines should copy and fill in).
- The `hPLZF` transgene contig is always appended to the mm10 + GENCODE vM25
  reference — any pipeline that quantifies gene expression for this project must
  include it, or the transgene signal (the primary biological readout) is lost.

## Environments

- Every Nextflow process gets its own micromamba environment YAML under the
  fork's `envs/` directory — no shared global environment, no bare `conda install`
  at the top level. Pin exact versions; for `forks/star_featurecounts/`, versions
  must match Plasmidsaurus's `BM9FTP-methods.txt` exactly (that's the point of the
  fork — a faithful independent reproduction).

## Architecture decisions

- Any non-obvious methodological or architectural choice gets an ADR in
  `docs/adr/` (use `docs/adr/template.md`). Examples already on file: why
  Nextflow, why micromamba, why the hybrid reference, why UMI dedup only applies
  to the STAR fork and not salmon/kallisto.
- Update `logs/session-log.md` at the end of each work session: what was done,
  why, and what's left.

## Adding a new analytical fork

This project is expected to grow multiple analytical forks that ask related but
independent questions (not necessarily sharing the same data or pipeline). Use the
`new-analytical-fork` skill (`.claude/skills/new-analytical-fork/SKILL.md`) to
scaffold a new fork directory consistent with the existing `star_featurecounts` /
`salmon_kallisto` pattern, rather than improvising a new structure per fork.
