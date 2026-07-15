// RSeQC v5.0.4 infer_experiment.py (part of methods.txt step 6): determines
// library strandedness per sample from the deduplicated BAM + the combined
// BED12 gene model, rather than hardcoding a strand assumption. The parsed
// result (a featureCounts -s code: 0 = unstranded, 1 = stranded, 2 = reverse
// stranded) is what main.nf uses -- by majority vote across all samples --
// to set the strand flag for the single batched featureCounts run.

process INFER_STRANDEDNESS {
    tag "${sample_id}"

    conda "${moduleDir}/../envs/rseqc.yml"

    publishDir { "${params.outdir}/rseqc/${sample_id}" }, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path bed12

    output:
    tuple val(sample_id), path("${sample_id}.infer_experiment.txt"), emit: report
    tuple val(sample_id), path("${sample_id}.strand_code.txt"),      emit: strand_code

    script:
    """
    infer_experiment.py -i ${bam} -r ${bed12} > ${sample_id}.infer_experiment.txt

    python3 - "${sample_id}.infer_experiment.txt" "${sample_id}.strand_code.txt" "${params.strand_unstranded_tolerance}" <<'PYEOF'
import sys

report_path, out_path, tolerance = sys.argv[1], sys.argv[2], float(sys.argv[3])

frac_forward = frac_reverse = None
with open(report_path) as fh:
    for line in fh:
        if '"++,--"' in line:
            frac_forward = float(line.strip().split(":")[-1])
        elif '"+-,-+"' in line:
            frac_reverse = float(line.strip().split(":")[-1])

if frac_forward is None or frac_reverse is None:
    code = 0
elif abs(frac_forward - frac_reverse) < tolerance:
    code = 0
elif frac_forward > frac_reverse:
    code = 1
else:
    code = 2

with open(out_path, "w") as out:
    out.write(f"{code}\\n")
PYEOF
    """

    stub:
    """
    touch ${sample_id}.infer_experiment.txt
    echo 0 > ${sample_id}.strand_code.txt
    """
}
