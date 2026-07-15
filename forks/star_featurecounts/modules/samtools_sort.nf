// samtools v1.22.1 coordinate sort of the STAR output BAM (methods.txt step 4),
// plus indexing so downstream RSeQC/Qualimap/UMICollapse steps can random-access it.

process SAMTOOLS_SORT {
    tag "${sample_id}"

    conda "${moduleDir}/../envs/samtools.yml"

    publishDir { "${params.outdir}/star/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: sorted

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam ${bam}
    samtools index ${sample_id}.sorted.bam
    """

    stub:
    """
    touch ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai
    """
}
