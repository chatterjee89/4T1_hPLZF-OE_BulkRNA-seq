// UMICollapse v1.1.0 (methods.txt step 5): PCR/optical duplicate removal
// using the UMI parsed from the read-name suffix (e.g.
// "...:1833:1014_TAGGCGAACACCAA" -> UMI "TAGGCGAACACCAA"), NOT a separate
// UMI fastq/index read. UMICollapse's default `--umi-sep` is "_", which
// matches this read-name layout exactly, so it is passed explicitly below
// for clarity rather than relied on implicitly.
//
// Single-end BAM/SAM mode ("bam" command, no --paired) operates directly on
// the coordinate-sorted BAM from samtools_sort.nf.

process UMI_DEDUP {
    tag "${sample_id}"
    label 'big_mem'

    conda "${moduleDir}/../envs/umicollapse.yml"

    publishDir { "${params.outdir}/umi_dedup/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.bam.bai"), emit: dedup

    script:
    """
    umicollapse bam \\
      -i ${bam} \\
      -o ${sample_id}.dedup.bam \\
      -k ${params.umi_edit_distance} \\
      --umi-sep _

    samtools index ${sample_id}.dedup.bam
    """

    stub:
    """
    touch ${sample_id}.dedup.bam ${sample_id}.dedup.bam.bai
    """
}
