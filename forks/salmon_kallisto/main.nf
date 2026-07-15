#!/usr/bin/env nextflow
//
// Fork 2: salmon / kallisto -- orthogonal, alignment-free pseudoalignment
// quantification of the BM9FTP 4T1 hPLZF-OE bulk RNA-seq samples.
//
// Run from the repo root, e.g.:
//   nextflow run forks/salmon_kallisto/main.nf -profile micromamba
//
// See forks/salmon_kallisto/README.md and docs/adr/0004, 0005, 0006 for the
// reasoning behind reference construction, why this fork exists alongside
// forks/star_featurecounts/, and why no UMI dedup happens here.

nextflow.enable.dsl = 2

include { BUILD_TRANSCRIPTOME } from './modules/build_transcriptome.nf'
include { FASTP               } from './modules/fastp.nf'
include { SALMON_INDEX        } from './modules/salmon_index.nf'
include { SALMON_QUANT        } from './modules/salmon_quant.nf'
include { KALLISTO_INDEX      } from './modules/kallisto_index.nf'
include { KALLISTO_QUANT      } from './modules/kallisto_quant.nf'
include { TXIMPORT_TO_GENE    } from './modules/tximport_to_gene.nf'

workflow {

    // params.fastq_dir / mm10_fasta / mm10_gtf / hplzf_fasta come from
    // conf/local.config (gitignored; see conf/local.config.example). Fail
    // fast with a clear message rather than a confusing null-path error deep
    // in a process if that file hasn't been set up on this machine.
    def required_params = ['fastq_dir', 'mm10_fasta', 'mm10_gtf', 'hplzf_fasta', 'samplesheet']
    def missing = required_params.findAll { !params[it] }
    if (missing) {
        error "Missing required params: ${missing.join(', ')}. Copy conf/local.config.example to conf/local.config and fill in paths for this machine."
    }

    mm10_fasta  = file(params.mm10_fasta,  checkIfExists: true)
    mm10_gtf    = file(params.mm10_gtf,    checkIfExists: true)
    hplzf_fasta = file(params.hplzf_fasta, checkIfExists: true)

    BUILD_TRANSCRIPTOME(mm10_fasta, mm10_gtf, hplzf_fasta)

    // samplesheet's fastq_path column has a literal, unresolved $DATA_ROOT
    // placeholder (it documents where the file lives on the machine that
    // generated the sheet) -- we deliberately ignore it and instead build the
    // real path from params.fastq_dir (this machine's actual data root) +
    // sample_id + .fastq.gz, since every sample_id already matches its fastq
    // filename stem exactly.
    samples_ch = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def fastq = file("${params.fastq_dir}/${row.sample_id}.fastq.gz")
            tuple(row.sample_id, fastq)
        }

    FASTP(samples_ch)

    SALMON_INDEX(BUILD_TRANSCRIPTOME.out.transcriptome_fasta)
    KALLISTO_INDEX(BUILD_TRANSCRIPTOME.out.transcriptome_fasta)

    // .first() turns the single-item index channel into a value channel so it
    // can be reused across all 24 per-sample SALMON_QUANT / KALLISTO_QUANT
    // calls instead of being drained after the first one.
    SALMON_QUANT(FASTP.out.trimmed, SALMON_INDEX.out.index.first())
    KALLISTO_QUANT(FASTP.out.trimmed, KALLISTO_INDEX.out.index.first())

    salmon_dirs_ch   = SALMON_QUANT.out.quant_dir.map { sample_id, dir -> dir }.collect()
    kallisto_dirs_ch = KALLISTO_QUANT.out.quant_dir.map { sample_id, dir -> dir }.collect()

    TXIMPORT_TO_GENE(
        BUILD_TRANSCRIPTOME.out.tx2gene,
        salmon_dirs_ch,
        kallisto_dirs_ch
    )
}
