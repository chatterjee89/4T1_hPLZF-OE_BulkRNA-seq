// Summarizes both salmon's and kallisto's transcript-level quantifications to
// gene-level counts via Bioconductor's tximport, using the shared tx2gene map
// from BUILD_TRANSCRIPTOME. Emits two separate gene-count matrices (one per
// quantifier) rather than merging them -- keeping them separate is deliberate,
// so that "does salmon agree with kallisto" and "does either agree with
// Plasmidsaurus/the STAR fork" stay answerable as distinct questions downstream
// in analysis/.
//
// The actual work happens in bin/tximport_to_gene.R (auto-added to PATH by
// Nextflow because it lives in this pipeline's bin/ directory); this process
// just stages inputs into the task's working directory under names the script
// expects (tx2gene.tsv, <sample_id>_salmon/, <sample_id>_kallisto/).

process TXIMPORT_TO_GENE {
    tag 'tximport_to_gene'

    conda "${moduleDir}/../envs/tximport.yml"

    publishDir "${params.outdir}/salmon_kallisto/gene_counts", mode: 'copy'

    input:
    path tx2gene
    path 'salmon_quant_dirs/*'
    path 'kallisto_quant_dirs/*'

    output:
    path 'salmon_gene_counts.tsv',   emit: salmon_gene_counts
    path 'kallisto_gene_counts.tsv', emit: kallisto_gene_counts

    script:
    """
    # Nextflow stages each collected quant dir under salmon_quant_dirs/ and
    # kallisto_quant_dirs/ to keep the two tools' outputs from colliding on
    # disk; symlink their contents back up to cwd (as <sample_id>_salmon/,
    # <sample_id>_kallisto/) so tximport_to_gene.R can glob for them directly.
    for d in salmon_quant_dirs/*_salmon; do ln -s "\$d" "\$(basename "\$d")"; done
    for d in kallisto_quant_dirs/*_kallisto; do ln -s "\$d" "\$(basename "\$d")"; done

    tximport_to_gene.R
    """

    stub:
    """
    printf 'gene_id\\tgene_name\\thPLZF_stub\\n' > salmon_gene_counts.tsv
    printf 'hPLZF\\thPLZF\\t10\\n' >> salmon_gene_counts.tsv
    printf 'gene_id\\tgene_name\\thPLZF_stub\\n' > kallisto_gene_counts.tsv
    printf 'hPLZF\\thPLZF\\t10\\n' >> kallisto_gene_counts.tsv
    """
}
