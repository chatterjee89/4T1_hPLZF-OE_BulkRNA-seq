# 0007. edgeR + Hallmark GSEA (in R) for cross-method comparison

Date: 2026-07-14
Status: Accepted

## Context

To compare biological conclusions (not just raw counts) across the Plasmidsaurus,
Fork 1, and Fork 2 count matrices, we need a consistent differential expression and
pathway enrichment method applied identically to all three. Plasmidsaurus's own
methods used edgeR (TMM normalization, `filterByExpr` default filtering) for DE and
gseapy (Python) with MSigDB Hallmark gene sets for GSEA.

## Decision

- Differential expression: edgeR, TMM normalization, `filterByExpr` defaults,
  EV vs PLZF contrast fit per timepoint — matches Plasmidsaurus's method exactly
  for direct comparability.
- Pathway enrichment: Hallmark GSEA via `fgsea` in R rather than Python's gseapy.
  The gene sets and method (pre-ranked GSEA against MSigDB Hallmark) are
  equivalent; keeping it in R avoids introducing a second language runtime into
  the `analysis/` stage purely for one script, since edgeR/PCA/plotting are
  already R.

## Consequences

Analysis results should be comparable to Plasmidsaurus's own DE calls when run on
their matrix (a useful internal check that our re-implementation is faithful). If
`fgsea` output ever meaningfully diverges from what `gseapy` would produce (e.g.
differing GSEA implementation gsea details), revisit and note discrepancies in a
follow-up ADR rather than silently reconciling them.
