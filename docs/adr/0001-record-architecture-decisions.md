# 0001. Record architecture decisions

Date: 2026-07-14
Status: Accepted

## Context

This project involves a series of non-obvious, hard-to-reverse methodological
choices (reference construction, UMI handling, quantification method, statistical
methods) made under real constraints (matching an existing ground-truth analysis,
working with a mouse+human hybrid system). Future contributors (including future
us) need to know *why* a choice was made, not just what it is.

## Decision

Record every non-obvious architectural or methodological decision as an ADR in
`docs/adr/`, numbered sequentially, using `template.md`.

## Consequences

Slightly more overhead per decision. In exchange, the rationale behind e.g. "why
does the salmon/kallisto fork skip UMI dedup" survives past the session that made
the call.
