# Downstream analysis: sample characterization, DGE, GSEA

Runs identically against any of the three gene-level count-matrix sources
(Plasmidsaurus / Fork 1 STAR+featureCounts / Fork 2 salmon+kallisto), so biological
conclusions — not just raw counts — can be compared across quantification methods.

- `01_sample_characterization.R` — TMM-normalized PCA + sample-sample correlation
  heatmap, colored by condition (EV/PLZF) × timepoint (D0–D3); a QC gate before
  trusting any DGE downstream.
- `02_dge_edgeR.R` — edgeR (TMM norm, `filterByExpr` defaults), EV vs PLZF contrast
  per timepoint, mirroring Plasmidsaurus's method for direct comparability. Flags
  the `hPLZF` transgene row explicitly as a built-in sanity check (should be ~0 in
  EV, high in PLZF).
- `03_gsea_hallmark.R` — Hallmark GSEA (`fgsea`) on each timepoint's DGE ranks.
- `04_cross_method_concordance.R` — once Fork 1/Fork 2 counts exist, compares DGE
  hit lists and GSEA pathway hits across all three sources (Layer 2 robustness
  check — see README.md at repo root and `docs/adr/0007`).

Currently run only against the real Plasmidsaurus matrix (`conf/local.config` →
`params.plasmidsaurus_matrix`), since Fork 1/Fork 2 pipelines are built but not yet
executed. Output plots land in `plots/` (gitignored, regenerated on each run).

Status: scaffolding in progress (`R/utils.R` and numbered scripts to follow).
