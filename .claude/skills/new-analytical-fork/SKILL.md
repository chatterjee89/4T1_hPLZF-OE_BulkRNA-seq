---
name: new-analytical-fork
description: Scaffold a new analytical fork directory under forks/, matching the structure used by star_featurecounts and salmon_kallisto. Use when starting a new pipeline/analysis fork that asks a related but independent question and doesn't necessarily share data with existing forks.
---

# New analytical fork

This project (4T1_hPLZF-OE_BulkRNA-seq) is expected to grow multiple analytical
forks over time — pipelines or analyses that ask related but independent
questions, not necessarily sharing the same data or method as existing forks
(`forks/star_featurecounts/`, `forks/salmon_kallisto/`).

## When invoked

1. Ask the user (if not already clear from context): fork name (kebab_case or
   snake_case, matching existing `star_featurecounts` / `salmon_kallisto` style),
   one-sentence purpose, and whether it's a Nextflow pipeline (needs `modules/` +
   `envs/`) or a pure analysis fork (more like `analysis/` — R/Python scripts, no
   Nextflow needed).

2. Create `forks/<fork_name>/` with:
   - `README.md` — purpose, method/tools, how it relates to (or differs from)
     existing forks, and its status (scaffolding vs. implemented vs. executed).
   - If a Nextflow pipeline: `main.nf`, `modules/` (one `.nf` file per logical
     step), `envs/` (one pinned micromamba YAML per process — no shared/global
     env, per `CLAUDE.md`).
   - If a pure analysis fork: numbered scripts (`01_...`, `02_...`) following the
     convention in `analysis/`, plus a `plots/` dir (gitignored — add a
     `.gitkeep`).

3. If the fork introduces a non-obvious methodological choice (a new reference,
   a new statistical method, a divergence from how existing forks handle
   something like UMIs or normalization), write an ADR in `docs/adr/` using
   `docs/adr/template.md` — don't skip this even for "obvious in the moment"
   choices; see `docs/adr/0001` for why.

4. Add an entry to `logs/session-log.md` under today's date describing what was
   scaffolded and why.

5. Do not wire the new fork into any other fork's code or config. Forks are
   independent — a new fork should be addable/removable without touching
   `forks/star_featurecounts/` or `forks/salmon_kallisto/`.
