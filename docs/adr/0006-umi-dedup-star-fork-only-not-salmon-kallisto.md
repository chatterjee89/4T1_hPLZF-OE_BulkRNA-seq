# 0006. UMI deduplication applies only to the STAR fork, not salmon/kallisto

Date: 2026-07-14
Status: Accepted

## Context

Reads carry UMIs embedded in the read-name suffix (e.g.
`...:1833:1014_TAGGCGAACACCAA`). Plasmidsaurus deduplicated PCR/optical duplicates
using UMIcollapse *after* STAR alignment, before featureCounts. For alignment-free
pseudoalignment (salmon/kallisto), standard practice diverges: these tools already
model multi-mapping and fragment-level ambiguity probabilistically via an EM
algorithm, and removing UMI duplicates prior to pseudoalignment is not standard
practice for bulk RNA-seq — doing so can bias abundance estimates rather than
correct them, unlike the alignment-based case where duplicate removal is a
well-established step before simple read counting.

## Decision

- `forks/star_featurecounts/`: dedup with UMIcollapse (same tool/version as
  Plasmidsaurus) after alignment, before featureCounts — matches ground truth
  methodology exactly.
- `forks/salmon_kallisto/`: run salmon/kallisto directly on fastp-trimmed reads,
  with **no UMI-based deduplication** step.

## Consequences

This is a known, deliberate divergence between the two forks' inputs. When
comparing Fork 2's counts to Plasmidsaurus/Fork 1, some discordance driven purely
by this methodological difference (not a bug) should be expected, particularly for
highly-expressed genes where PCR duplication rates are higher. `analysis/
04_cross_method_concordance.R` should call this out explicitly rather than
treating all discordance as pipeline error.
