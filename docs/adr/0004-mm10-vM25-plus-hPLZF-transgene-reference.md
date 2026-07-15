# 0004. mm10 + GENCODE vM25 reference with an hPLZF transgene contig

Date: 2026-07-14
Status: Accepted

## Context

ZBTB16/PLZF is overexpressed in these 4T1 (mouse) cells via a *human* ZBTB16
(hPLZF) transgene construct. A standard mouse reference alone would either drop
reads originating from the transgene or, worse, silently misassign them to the
endogenous mouse `Zbtb16` locus (the human and mouse coding sequences are similar
enough to risk cross-mapping), corrupting both the transgene-expression readout and
the endogenous `Zbtb16` counts. Plasmidsaurus's existing expression matrix already
contains a `256578 hPLZF transgene` row, confirming they solved this the same way.

## Decision

Build the reference genome/transcriptome from UCSC mm10 (`mm10_main.fa`, 21 main
contigs, no random/alt scaffolds — matching what's already on disk) + GENCODE vM25
annotation, with the hPLZF cDNA (`reference/hZBTB16_cDNA/hPLZF.fa`) appended as an
extra contig/transcript and a corresponding synthetic GTF/gene entry, in both the
STAR fork (genome-level contig) and the salmon/kallisto fork (transcript-level
entry + `tx2gene` mapping).

## Consequences

Both forks must include a `build_reference` / `build_transcriptome` step rather
than pointing at an off-the-shelf reference. This keeps our results comparable to
Plasmidsaurus's (same hybrid reference strategy) and preserves the ability to
directly read off transgene expression instead of losing that signal or
contaminating the endogenous `Zbtb16` counts.
