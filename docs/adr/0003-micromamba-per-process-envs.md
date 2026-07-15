# 0003. micromamba per-process environments, not Docker

Date: 2026-07-14
Status: Accepted

## Context

Nothing was installed on this machine at project start (no Docker, no conda/mamba).
Each Nextflow process needs an isolated, reproducible environment (STAR, fastp,
samtools, UMIcollapse, subread, salmon, kallisto, RSeQC, Qualimap, MultiQC, R/
Bioconductor packages) without those tools' dependencies colliding.

## Decision

Use micromamba (a small, fast, daemon-less conda-compatible environment manager)
with one YAML environment spec per Nextflow process (`envs/*.yml` in each fork),
selected via Nextflow's `conda` profile — rather than Docker Desktop containers.

## Consequences

Lighter installation footprint (no background daemon, no ~700MB+ Docker Desktop
install, no license considerations) and standard practice in academic genomics.
Trade-off: slightly less portable/reproducible than fully pinned containers (conda
solves can drift over time even with pinned versions; no container-level OS
isolation). If this project is ever run on shared HPC infrastructure or needs
byte-for-byte reproducible containers, revisit this decision — Nextflow's `conda`
and `container` directives can coexist via profiles without restructuring the
pipelines.
