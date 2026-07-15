// salmon quant for single-end reads.
//
// ASSUMPTION TO REVISIT: these are single-end reads, so salmon cannot observe
// the true fragment-length distribution from read pairs. --fldMean 200/--fldSD 80
// are reasonable generic defaults for ~95bp SE Illumina RNA-seq libraries but are
// NOT measured from this data -- once real reads are run, sanity-check against
// e.g. RSeQC/Qualimap insert-size-adjacent metrics from the STAR fork's BAMs (or
// vendor library-prep documentation, if available) and update here if warranted.
// kallisto_quant.nf uses the identical assumption so the two tools stay comparable.

process SALMON_QUANT {
    tag "$sample_id"
    // Pseudoalignment quant against a transcriptome index is lighter than STAR
    // genome alignment -- default base.config memory plus a modest cpu bump for
    // salmon's threaded quantification, not the shared `alignment` label.
    cpus 4

    conda "${moduleDir}/../envs/salmon.yml"

    publishDir "${params.outdir}/salmon_kallisto/salmon/quant", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_fastq)
    path salmon_index

    output:
    tuple val(sample_id), path("${sample_id}_salmon"), emit: quant_dir

    script:
    """
    salmon quant \\
        -i ${salmon_index} \\
        -l A \\
        -r ${trimmed_fastq} \\
        --validateMappings \\
        --fldMean 200 \\
        --fldSD 80 \\
        -p ${task.cpus} \\
        -o ${sample_id}_salmon
    """

    stub:
    """
    mkdir -p ${sample_id}_salmon
    printf 'Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\n' > ${sample_id}_salmon/quant.sf
    printf 'hPLZF\\t1000\\t900\\t1.0\\t10\\n' >> ${sample_id}_salmon/quant.sf
    """
}
