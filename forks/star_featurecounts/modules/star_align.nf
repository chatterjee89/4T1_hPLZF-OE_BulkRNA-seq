// STAR v2.7.11(b) alignment of fastp-trimmed reads against the combined
// mm10+hPLZF index. Non-canonical splice junctions are removed
// (--outFilterIntronMotifs RemoveNoncanonical) and unmapped reads are written
// out as fastx (--outReadsUnmapped Fastx), matching methods.txt step 3.
// Coordinate sorting is deliberately NOT done here -- that's its own step
// (samtools_sort.nf), mirroring the vendor's method as a distinct step 4.

process STAR_ALIGN {
    tag "${sample_id}"
    label 'alignment'

    conda "${moduleDir}/../envs/star_align.yml"

    publishDir { "${params.outdir}/star/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(fastq)
    path star_index

    output:
    tuple val(sample_id), path("${sample_id}.Aligned.out.bam"),        emit: bam
    tuple val(sample_id), path("${sample_id}.Unmapped.out.mate1"),     emit: unmapped, optional: true
    tuple val(sample_id), path("${sample_id}.Log.final.out"),          emit: log_final
    tuple val(sample_id), path("${sample_id}.SJ.out.tab"),             emit: splice_junctions

    script:
    """
    STAR \\
      --runMode alignReads \\
      --genomeDir ${star_index} \\
      --readFilesIn ${fastq} \\
      --readFilesCommand zcat \\
      --outFilterIntronMotifs RemoveNoncanonical \\
      --outReadsUnmapped Fastx \\
      --outSAMtype BAM Unsorted \\
      --outFileNamePrefix ${sample_id}. \\
      --runThreadN ${task.cpus}
    """

    stub:
    """
    touch ${sample_id}.Aligned.out.bam
    touch ${sample_id}.Unmapped.out.mate1
    touch ${sample_id}.Log.final.out
    touch ${sample_id}.SJ.out.tab
    """
}
