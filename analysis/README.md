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
executed. Output plots land in `plots/`; DGE and GSEA result tables land in
`results/` (both regenerated on each run — see "Output locations" below).

Status: all four scripts are implemented and have been run end-to-end against
the real Plasmidsaurus matrix. `R/utils.R` provides the shared loading/joining/
theming helpers; `04_cross_method_concordance.R` runs in single-source mode
until Fork 1/Fork 2 produce real output (see below).

## Layout

- `R/utils.R` — shared helpers:
  - `load_count_matrix(path, source)` — loads a count matrix into a common
    long-form tibble (`gene_id, gene_name, gene_biotype, sample_num, count`).
    Only `source = "plasmidsaurus"` is implemented today; `"star_featurecounts"`
    and `"salmon_kallisto"` are stubbed with a clear error (not fabricated data)
    describing the expected output shape, to be filled in once those forks run.
  - `read_sample_sheet()` — reads `samplesheet/BM9FTP_samplesheet.csv`.
  - `build_count_matrix()` — joins counts to sample sheet metadata **via
    `sample_num`, not `sample_id`** — the Plasmidsaurus matrix's sample columns
    are labeled `BM9FTP_<sample_num>_{cpm,count}`, not the full `sample_id`
    (which also carries condition/timepoint/rep), so a direct string join
    would silently fail to match. Also builds a `genes x samples` numeric
    matrix plus a matching per-sample metadata table.
  - `theme_analysis()`, `condition_colors` — shared ggplot2 look, consistent
    across all plots.
  - `analysis_plots_dir()` / `analysis_results_dir()` — resolve/create
    `analysis/plots/` and `analysis/results/` (scripts are run from repo root).
  - `HPLZF_GENE_ID <- "256578"` — single source of truth for the transgene
    sanity-check row, used by scripts 01/02/04.
- `01_sample_characterization.R` — DGEList + TMM normalization, PCA on
  log-CPM (scatter colored by condition, faceted by timepoint), and a
  sample-sample Pearson correlation heatmap (`pheatmap`, annotated by
  condition + timepoint).
- `02_dge_edgeR.R` — edgeR DE, **fit separately per timepoint** (4 independent
  3-vs-3 EV-vs-PLZF contrasts: D0, D1, D2, D3), TMM normalization,
  `filterByExpr` defaults, `glmQLFit`/`glmQLFTest` (see "Design decisions").
  One result table + one volcano plot per timepoint; explicit hPLZF
  sanity-check summary printed and written to
  `results/02_hplzf_sanity_check_<source>.tsv`.
- `03_gsea_hallmark.R` — pre-ranked GSEA (`fgsea`) against MSigDB Hallmark
  gene sets (`msigdbr`, mouse-native), ranked per-timepoint on each DGE
  result (see "Design decisions" for the ranking metric). One result table +
  one dotplot (top 15 pathways by padj) per timepoint.
- `04_cross_method_concordance.R` — designed to compare DGE hit-list overlap
  (Jaccard) and GSEA pathway hits across however many of
  {plasmidsaurus, star_featurecounts, salmon_kallisto} have real
  `02_dge_edgeR.R` output on disk. Today only `plasmidsaurus` qualifies, so it
  runs in **single-source mode**: it summarizes/validates that one source
  (gene counts tested, # significant, hPLZF logFC/FDR per timepoint) and
  prints an explicit notice that Fork 1/Fork 2 comparison will activate
  automatically, with no fabricated placeholder data, once those pipelines
  produce real output.

## Output locations

- `plots/` — all PNGs from scripts 01–03 (gitignored; regenerated per run).
- `results/` — TSV result tables from scripts 02–04 (DGE tables, hPLZF
  sanity-check summary, GSEA tables, single/multi-source concordance
  summaries). **This directory is new** (not present in the original stub)
  and, like `plots/`, is regenerated output rather than source — it will need
  a `.gitignore` entry (e.g. `analysis/results/*` with a `.gitkeep`
  exception, mirroring the existing `analysis/plots/*` pattern) added at the
  repo root; this analysis agent did not touch `.gitignore` since it lives
  outside `analysis/`.

All scripts accept an optional `--source=<name>` argument (default
`plasmidsaurus`); file names are suffixed with the source so multiple
sources' outputs can coexist once Fork 1/Fork 2 are run.

## Design decisions

- **DGE test: `glmQLFit`/`glmQLFTest`, not `exactTest`.** Quasi-likelihood
  F-tests give more robust error control than `exactTest` at these small
  (n=3 per group per timepoint) sample sizes, via empirical Bayes moderation
  of the QL dispersion, and this is the modern recommended edgeR workflow.
  `filterByExpr` is run **separately within each timepoint's 3-vs-3 subset**
  (not once genome-wide), which matters for the hPLZF sanity check — see
  below.
- **GSEA ranking metric: `-log10(PValue) * sign(logFC)`.** Chosen over plain
  `logFC` so genes are ranked by evidence strength *and* direction rather
  than by (sometimes noisy) fold-change magnitude alone. Documented here per
  `docs/adr/0007`; flagged as revisitable if exact parity with Plasmidsaurus's
  `gseapy` ranking is later required.
- **Hallmark gene sets via `msigdbr`.** This `msigdbr` installation (v26.1.0)
  uses a newer API than the classic `category = "H"` / `Mm.H` convention: mouse
  Hallmark sets are fetched via `msigdbr(db_species = "MM", species = "Mus musculus",
  collection = "MH")` (confirmed against `msigdbr_collections(db_species = "MM")`).
  Functionally this is the same mouse-native Hallmark collection (50 gene sets)
  the older `Mm.H` naming referred to.

## hPLZF transgene sanity check — result

The whole point of wiring `sample_num`-based joins and per-timepoint filtering
correctly is that this number should come out right. It did:

| Timepoint | logFC (PLZF vs EV) | PValue | FDR |
|---|---|---|---|
| D0 | *filtered out* (near-zero counts in both arms; see below) | — | — |
| D1 | 6.92 | 2.1e-08 | 3.5e-06 |
| D2 | 7.97 | 5.3e-10 | 2.1e-07 |
| D3 | 4.61 | 5.7e-07 | 8.5e-05 |

Raw hPLZF counts confirm this isn't a filtering artifact: EV samples read 0–2
counts at every timepoint (background/noise level), while PLZF samples read
~0 at D0 (transgene not yet induced) then jump to 11–40 counts from D1
onward. At D0, `filterByExpr` correctly drops the gene from *both* the
per-timepoint D0 contrast (only 1/3 PLZF replicates non-zero, "expression in
all 3 replicates" default threshold not met — i.e., biologically absent, not
a bug) and, separately, from the genome-wide QC pass in
`01_sample_characterization.R` (there, with 12-vs-12 samples exceeding
`filterByExpr`'s `large.n = 10`, the default `min.prop = 0.7` inflates the
required expressing-sample count to ~12/12, so a gene induced only from D1
onward fails that pass too — expected given the biology, not a wiring bug).
From D1 onward hPLZF is the strongest logFC hit in every volcano plot (see
`plots/02_volcano_D1..D3_plasmidsaurus.png`, gold-highlighted point), exactly
the induction-kinetics pattern expected of a transgene that isn't fully
switched on until after D0. This is strong evidence the sample_num-based
join, per-timepoint contrast design, and DGE pipeline are all wired correctly.

GSEA at D3 (furthest timepoint, largest effect) shows the expected
proliferation/growth signature: `E2F_TARGETS`, `G2M_CHECKPOINT`, and
`MYC_TARGETS_V1` are the top up-in-PLZF Hallmark pathways (all padj < 0.002),
with `EPITHELIAL_MESENCHYMAL_TRANSITION` and `HYPOXIA` among the top
down-in-PLZF pathways — directionally consistent with a proliferative,
less-hypoxic/EMT phenotype under PLZF overexpression, though a full
biological interpretation is out of scope for this scaffolding pass.

## Possible follow-up (flagged, not written here)

Whether the genome-wide QC-pass `filterByExpr` behavior in
`01_sample_characterization.R` (dropping hPLZF due to the `large.n`/`min.prop`
interaction, distinct from the per-timepoint contrasts in `02_dge_edgeR.R`
where it behaves as expected) merits its own ADR, or is adequately covered by
the inline comments/warnings already in the script, is left for the
orchestrating session to decide.
