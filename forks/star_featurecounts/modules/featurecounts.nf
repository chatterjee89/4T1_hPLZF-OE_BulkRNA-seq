// featureCounts (subread v2.1.1), methods.txt step 8: strand-specific
// counting, multi-mapping reads counted with fractional assignment
// (-M --fraction), feature types exon AND three_prime_UTR (-t
// exon,three_prime_UTR -g gene_id).
//
// Run once, batched across all 24 deduplicated BAMs (see fork README for why
// batch vs. per-sample was chosen) so the output is a single gene x sample
// count matrix directly comparable in shape to the Plasmidsaurus ground-truth
// matrix. Column headers in the output are the staged BAM filenames
// ("<sample_id>.dedup.bam"); downstream code strips the ".dedup.bam" suffix
// to recover sample_id.
//
// Strand code (0/1/2) comes from main.nf's majority vote over all
// per-sample infer_experiment.py results (infer_strandedness.nf) -- inferred
// from the data rather than hardcoded, per the fork's design brief.

process FEATURECOUNTS {
    tag "all_samples"
    label 'big_mem'

    conda "${moduleDir}/../envs/featurecounts.yml"

    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    path bams
    path gtf
    val strand_code

    output:
    path 'gene_counts.tsv',         emit: counts
    path 'gene_counts.tsv.summary', emit: summary

    script:
    """
    featureCounts \\
      -a ${gtf} \\
      -o gene_counts.tsv \\
      -t exon,three_prime_UTR \\
      -g gene_id \\
      -M --fraction \\
      -s ${strand_code} \\
      -T ${task.cpus} \\
      ${bams}
    """

    stub:
    """
    touch gene_counts.tsv gene_counts.tsv.summary
    """
}
