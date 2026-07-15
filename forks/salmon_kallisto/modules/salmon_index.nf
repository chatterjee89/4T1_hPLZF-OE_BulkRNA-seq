// Builds a salmon index from the mm10 + hPLZF transcript FASTA.
//
// NOTE: this is a plain (non-decoy-aware) transcriptome index -- it does not
// include the mm10 genome as a decoy sequence. Decoy-aware indexing is salmon's
// recommended best practice (reduces spurious mapping of intergenic/intronic
// fragments to transcripts with partial sequence similarity) but requires
// indexing against the full ~2.7GB mm10 genome as well as the transcriptome,
// which is a much larger one-time cost. Flagged here as a documented
// simplification / candidate follow-up once this fork is actually run.

process SALMON_INDEX {
    tag 'salmon_index'
    // Indexing a transcriptome (not a genome) is far lighter than STAR genome
    // generation, so this deliberately does not use the shared `alignment`
    // label (32GB, sized for the sibling STAR fork) -- just bump cpus for
    // salmon's multithreaded indexing and keep base.config's default memory.
    cpus 4

    conda "${moduleDir}/../envs/salmon.yml"

    publishDir "${params.outdir}/salmon_kallisto/salmon", mode: 'copy'

    input:
    path transcriptome_fasta

    output:
    path 'salmon_index', emit: index

    script:
    """
    salmon index \\
        -t ${transcriptome_fasta} \\
        -i salmon_index \\
        -k 31 \\
        -p ${task.cpus}
    """

    stub:
    """
    mkdir -p salmon_index
    touch salmon_index/info.json
    """
}
