// Builds a kallisto index from the same mm10 + hPLZF transcript FASTA used for
// salmon (modules/build_transcriptome.nf), so both tools quantify against an
// identical set of transcript sequences.

process KALLISTO_INDEX {
    tag 'kallisto_index'
    // kallisto index over a transcriptome is single-threaded and lighter than
    // STAR genome generation -- default base.config resources are enough, no
    // need for the shared `alignment` label (32GB, sized for the STAR fork).

    conda "${moduleDir}/../envs/kallisto.yml"

    publishDir "${params.outdir}/salmon_kallisto/kallisto", mode: 'copy'

    input:
    path transcriptome_fasta

    output:
    path 'kallisto_index.idx', emit: index

    script:
    """
    kallisto index -i kallisto_index.idx ${transcriptome_fasta}
    """

    stub:
    """
    touch kallisto_index.idx
    """
}
