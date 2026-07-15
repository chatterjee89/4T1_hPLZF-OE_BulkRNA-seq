// MultiQC v1.32 (methods.txt step 7): aggregates FastQC, fastp, STAR,
// RSeQC (infer_experiment + read_distribution), Qualimap, and featureCounts
// reports into one HTML report.

process MULTIQC {
    tag "multiqc"

    conda "${moduleDir}/../envs/multiqc.yml"

    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path reports

    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_data',        emit: data

    script:
    """
    multiqc . --force --filename multiqc_report.html
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_data
    touch multiqc_data/placeholder.txt
    """
}
