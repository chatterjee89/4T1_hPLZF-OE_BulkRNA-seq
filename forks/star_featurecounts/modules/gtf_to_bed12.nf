// Converts the combined GTF (mm10 + GENCODE vM25 + hPLZF) to a BED12 file,
// which RSeQC's infer_experiment.py and read_distribution.py require as
// their gene-model input (they don't accept GTF directly). Run once against
// the combined annotation, not per-sample.
//
// Not one of Plasmidsaurus's pinned tools -- see envs/ucsc_tools.yml for why
// these UCSC "kent" utilities are unpinned relative to a vendor version.

process GTF_TO_BED12 {
    tag "combined_annotation"

    conda "${moduleDir}/../envs/ucsc_tools.yml"

    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path gtf

    output:
    path 'combined_annotation.bed12', emit: bed12

    script:
    """
    gtfToGenePred ${gtf} combined_annotation.genePred
    genePredToBed combined_annotation.genePred combined_annotation.bed12
    """

    stub:
    """
    touch combined_annotation.bed12
    """
}
