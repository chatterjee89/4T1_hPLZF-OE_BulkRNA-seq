// kallisto quant for single-end reads.
//
// Unlike salmon, kallisto REQUIRES --fragment-length/--sd for single-end input
// (it has no other way to estimate effective transcript length). Uses the same
// mean=200/sd=80 assumption as salmon_quant.nf, for the same "revisit once real
// data is run" reason documented there -- keeping both tools on one shared
// assumption also keeps their outputs comparable to each other, not just to
// ground truth.
//
// --plaintext requests abundance.tsv (instead of abundance.h5) so downstream
// tximport (modules/tximport_to_gene.nf) doesn't need an HDF5 R dependency.

process KALLISTO_QUANT {
    tag "$sample_id"
    // Same reasoning as salmon_quant.nf: pseudoalignment quant is lighter than
    // STAR alignment, so this uses a modest cpu bump rather than the shared
    // `alignment` label sized for the STAR fork.
    cpus 4

    conda "${moduleDir}/../envs/kallisto.yml"

    publishDir "${params.outdir}/salmon_kallisto/kallisto/quant", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_fastq)
    path kallisto_index

    output:
    tuple val(sample_id), path("${sample_id}_kallisto"), emit: quant_dir

    script:
    """
    kallisto quant \\
        -i ${kallisto_index} \\
        -o ${sample_id}_kallisto \\
        --single \\
        --fragment-length 200 \\
        --sd 80 \\
        --plaintext \\
        -t ${task.cpus} \\
        ${trimmed_fastq}
    """

    stub:
    """
    mkdir -p ${sample_id}_kallisto
    printf 'target_id\\tlength\\teff_length\\test_counts\\ttpm\\n' > ${sample_id}_kallisto/abundance.tsv
    printf 'hPLZF\\t1000\\t900\\t10\\t1.0\\n' >> ${sample_id}_kallisto/abundance.tsv
    touch ${sample_id}_kallisto/run_info.json
    """
}
