# 0002. Nextflow for pipeline orchestration

Date: 2026-07-14
Status: Accepted

## Context

This project will contain multiple analytical forks (STAR-based, salmon/kallisto-
based, and future ones) that share some steps (QC, trimming) but diverge in
alignment/quantification strategy, tool versions, and resource needs. We need an
orchestration layer that: manages per-step software environments cleanly, resumes
from partial runs (large genomics jobs fail partway through), and is a domain
standard so the pipelines are legible to other bioinformaticians.

## Decision

Use Nextflow (DSL2) to orchestrate each fork. Snakemake and plain scripts/Makefile
were considered; Nextflow was chosen for its strong per-process environment
isolation (`conda`/`container` directives scoped to individual processes, which
matters here since Fork 1 and Fork 2 have almost entirely disjoint toolchains) and
built-in resumability (`-resume`), which matters given multi-hour STAR alignment
steps over 24 large fastq files.

## Consequences

Adds a JVM dependency (Java) to the toolchain, and a learning curve for anyone not
already familiar with Nextflow/Groovy-ish DSL. In exchange, each fork's environment
declarations stay local to the fork rather than requiring one global environment
that satisfies every tool's dependency constraints simultaneously.
