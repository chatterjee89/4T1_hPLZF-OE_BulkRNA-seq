#!/usr/bin/env nextflow
/*
 * Fork 1: STAR + featureCounts
 *
 * Independently rebuilds Plasmidsaurus's quantification method from raw
 * fastq (not their BAMs), pinned to their exact tool versions, so the
 * resulting gene x sample count matrix is directly comparable to their
 * ground-truth expression matrix.
 *
 * Pipeline: build_reference -> fastqc (informational) -> fastp -> STAR ->
 * samtools sort -> UMICollapse dedup -> RSeQC/Qualimap QC + strandedness
 * inference -> featureCounts (batched) -> MultiQC.
 *
 * See docs/adr/0004 (hybrid reference), 0005 (why two forks), 0006 (UMI
 * dedup scope) for the reasoning behind the non-obvious choices here, and
 * this fork's README.md for tool-version substitutions and design notes.
 */

nextflow.enable.dsl = 2

include { BUILD_REFERENCE }         from './modules/build_reference.nf'
include { FASTQC_RAW }              from './modules/fastqc.nf'
include { FASTP }                   from './modules/fastp.nf'
include { STAR_ALIGN }              from './modules/star_align.nf'
include { SAMTOOLS_SORT }           from './modules/samtools_sort.nf'
include { UMI_DEDUP }               from './modules/umi_dedup.nf'
include { GTF_TO_BED12 }            from './modules/gtf_to_bed12.nf'
include { INFER_STRANDEDNESS }      from './modules/infer_strandedness.nf'
include { RSEQC_READ_DISTRIBUTION } from './modules/rseqc_read_distribution.nf'
include { QUALIMAP_RNASEQ }         from './modules/qualimap.nf'
include { FEATURECOUNTS }           from './modules/featurecounts.nf'
include { MULTIQC }                 from './modules/multiqc.nf'

// Fork-specific defaults. Kept local to this fork's main.nf (rather than in
// the shared top-level nextflow.config / conf/*.config) since they're
// star_featurecounts-specific method parameters, not machine-specific paths.
// Override on the CLI, e.g. --fastp_min_length 30, if ever needed.
params.star_sjdb_overhang          = params.star_sjdb_overhang          ?: 94   // read length (~95bp) - 1, per STAR manual
params.fastp_min_quality           = params.fastp_min_quality           ?: 15   // methods.txt: min Phred quality 15
params.fastp_min_length            = params.fastp_min_length            ?: 50   // methods.txt: min length 50bp
params.umi_edit_distance           = params.umi_edit_distance           ?: 1    // UMICollapse -k default
params.strand_unstranded_tolerance = params.strand_unstranded_tolerance ?: 0.15 // |frac_fwd - frac_rev| below this => unstranded

workflow {

    if( !params.mm10_fasta || !params.mm10_gtf || !params.hplzf_fasta || !params.fastq_dir ) {
        error "params.mm10_fasta / mm10_gtf / hplzf_fasta / fastq_dir must be set -- copy conf/local.config.example to conf/local.config and fill in paths for this machine."
    }

    // ---- 0. Hybrid reference (mm10 + GENCODE vM25 + hPLZF transgene) ----
    BUILD_REFERENCE(
        file(params.mm10_fasta),
        file(params.mm10_gtf),
        file(params.hplzf_fasta)
    )
    GTF_TO_BED12(BUILD_REFERENCE.out.gtf)

    // ---- 1. Samples, from the samplesheet (fastq path built from
    //         params.fastq_dir + sample_id, NOT from the samplesheet's
    //         literal "$DATA_ROOT" placeholder) ----
    samples_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def fastq = file("${params.fastq_dir}/${row.sample_id}.fastq.gz")
            tuple(row.sample_id, fastq)
        }

    // ---- 2. QC + trimming ----
    FASTQC_RAW(samples_ch)
    FASTP(samples_ch)

    // ---- 3. Alignment ----
    STAR_ALIGN(FASTP.out.trimmed, BUILD_REFERENCE.out.index)

    // ---- 4. Coordinate sort ----
    SAMTOOLS_SORT(STAR_ALIGN.out.bam)

    // ---- 5. UMI-based PCR/optical duplicate removal ----
    UMI_DEDUP(SAMTOOLS_SORT.out.sorted)

    // ---- 6. Strandedness inference + informational alignment QC, all
    //         against the deduplicated BAM ----
    INFER_STRANDEDNESS(UMI_DEDUP.out.dedup, GTF_TO_BED12.out.bed12)
    RSEQC_READ_DISTRIBUTION(UMI_DEDUP.out.dedup, GTF_TO_BED12.out.bed12)
    QUALIMAP_RNASEQ(UMI_DEDUP.out.dedup, BUILD_REFERENCE.out.gtf)

    // Majority vote across all 24 samples' inferred strand codes -> one
    // featureCounts -s value for the single batched counting run. (All
    // samples share the same library prep, so they're expected to agree;
    // the vote is a safety net against a stray borderline sample.)
    strand_code_ch = INFER_STRANDEDNESS.out.strand_code
        .map { sample_id, f -> f.text.trim().toInteger() }
        .collect()
        .map { codes ->
            def counts = [:]
            codes.each { c -> counts[c] = (counts[c] ?: 0) + 1 }
            counts.max { it.value }.key
        }

    // ---- 7. Counting, batched across all samples (see featurecounts.nf
    //         header comment for why batch rather than per-sample) ----
    dedup_bams_ch = UMI_DEDUP.out.dedup
        .map { sample_id, bam, bai -> bam }
        .collect()

    FEATURECOUNTS(dedup_bams_ch, BUILD_REFERENCE.out.gtf, strand_code_ch)

    // ---- 8. Aggregate QC ----
    multiqc_inputs_ch = FASTQC_RAW.out.zip.map { it[1] }
        .mix(FASTP.out.json.map { it[1] })
        .mix(STAR_ALIGN.out.log_final.map { it[1] })
        .mix(INFER_STRANDEDNESS.out.report.map { it[1] })
        .mix(RSEQC_READ_DISTRIBUTION.out.report.map { it[1] })
        .mix(QUALIMAP_RNASEQ.out.report.map { it[1] })
        .mix(FEATURECOUNTS.out.summary)
        .collect()

    MULTIQC(multiqc_inputs_ch)
}
