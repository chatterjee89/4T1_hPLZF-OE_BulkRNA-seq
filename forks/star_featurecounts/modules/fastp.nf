// fastp v0.24.0, single-end mode: poly-X tail trimming, 3' quality-based tail
// trimming, minimum Phred quality 15, minimum length 50bp (methods.txt step 2).

process FASTP {
    tag "${sample_id}"

    conda "${moduleDir}/../envs/fastp.yml"

    publishDir { "${params.outdir}/fastp/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: trimmed
    tuple val(sample_id), path("${sample_id}.fastp.json"),       emit: json
    tuple val(sample_id), path("${sample_id}.fastp.html"),       emit: html

    script:
    """
    fastp \\
      --in1 ${fastq} \\
      --out1 ${sample_id}.trimmed.fastq.gz \\
      --trim_poly_x \\
      --cut_tail \\
      --qualified_quality_phred ${params.fastp_min_quality} \\
      --length_required ${params.fastp_min_length} \\
      --thread ${task.cpus} \\
      --json ${sample_id}.fastp.json \\
      --html ${sample_id}.fastp.html
    """

    stub:
    """
    echo "" | gzip > ${sample_id}.trimmed.fastq.gz
    touch ${sample_id}.fastp.json ${sample_id}.fastp.html
    """
}
