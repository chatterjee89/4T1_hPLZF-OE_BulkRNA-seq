// RSeQC v5.0.4 read_distribution.py (methods.txt step 6): read distribution
// across genomic features (exons, introns, UTRs, intergenic). QC only,
// informational -- does not feed featureCounts, only MultiQC.

process RSEQC_READ_DISTRIBUTION {
    tag "${sample_id}"

    conda "${moduleDir}/../envs/rseqc.yml"

    publishDir { "${params.outdir}/rseqc/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path bed12

    output:
    tuple val(sample_id), path("${sample_id}.read_distribution.txt"), emit: report

    script:
    """
    read_distribution.py -i ${bam} -r ${bed12} > ${sample_id}.read_distribution.txt
    """

    stub:
    """
    touch ${sample_id}.read_distribution.txt
    """
}
