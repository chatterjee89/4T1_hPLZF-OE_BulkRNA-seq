// Trim-only QC on the raw single-end reads. Deliberately does NOT do any
// UMI-aware deduplication (the UMI lives as a suffix on the read name, e.g.
// ...:1833:1014_TAGGCGAACACCAA) -- see docs/adr/0006. salmon/kallisto's EM-based
// abundance estimation already models multi-mapping/fragment ambiguity, and
// pre-dedup on pseudoalignment input is not standard bulk RNA-seq practice.

process FASTP {
    tag "$sample_id"

    conda "${moduleDir}/../envs/fastp.yml"

    publishDir "${params.outdir}/salmon_kallisto/fastp", mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("${sample_id}.trimmed.fastq.gz"), emit: trimmed
    path "${sample_id}.fastp.json", emit: json
    path "${sample_id}.fastp.html", emit: html

    script:
    """
    fastp \\
        --in1 ${fastq} \\
        --out1 ${sample_id}.trimmed.fastq.gz \\
        --trim_poly_x \\
        --qualified_quality_phred 20 \\
        --length_required 36 \\
        --thread ${task.cpus} \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html
    """

    stub:
    """
    printf '' | gzip > ${sample_id}.trimmed.fastq.gz
    touch ${sample_id}.fastp.json ${sample_id}.fastp.html
    """
}
