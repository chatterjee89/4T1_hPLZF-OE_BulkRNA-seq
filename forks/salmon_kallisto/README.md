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

Status: scaffolding in progress (`main.nf` / `modules/` / `envs/` to follow).
