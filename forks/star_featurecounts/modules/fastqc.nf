// FastQC v0.12.1 on the raw fastq -- QC only, informational (methods.txt step
// 1). Doesn't gate or alter anything downstream; its report is folded into
// the final MultiQC report.

process FASTQC_RAW {
    tag "${sample_id}"

    conda "${moduleDir}/../envs/fastqc.yml"

    publishDir { "${params.outdir}/fastqc_raw/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("*_fastqc.zip"),  emit: zip
    tuple val(sample_id), path("*_fastqc.html"), emit: html

    script:
    """
    fastqc --threads ${task.cpus} ${fastq}
    """

    stub:
    """
    touch ${sample_id}_fastqc.zip ${sample_id}_fastqc.html
    """
}
