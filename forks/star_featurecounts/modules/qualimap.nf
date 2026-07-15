// Qualimap v2.3 rnaseq mode (methods.txt step 6): alignment/read-count QC
// against the combined annotation. QC only, informational -- does not feed
// featureCounts, only MultiQC. Single-end data, so no read-name sorting
// requirement (that only matters for Qualimap's paired-end pairing logic).

process QUALIMAP_RNASEQ {
    tag "${sample_id}"
    label 'big_mem'

    conda "${moduleDir}/../envs/qualimap.yml"

    publishDir { "${params.outdir}/qualimap/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path gtf

    output:
    tuple val(sample_id), path("${sample_id}_qualimap"), emit: report

    script:
    """
    qualimap rnaseq \\
      -bam ${bam} \\
      -gtf ${gtf} \\
      -outdir ${sample_id}_qualimap \\
      --java-mem-size=${task.memory.toGiga()}G
    """

    stub:
    """
    mkdir -p ${sample_id}_qualimap
    touch ${sample_id}_qualimap/rnaseq_qc_results.txt
    """
}
