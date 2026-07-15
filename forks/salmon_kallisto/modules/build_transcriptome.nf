// Builds the transcript-level reference used by both salmon and kallisto:
//   1. Extract spliced transcript FASTA from mm10 (main contigs only) + GENCODE
//      vM25 GTF via gffread.
//   2. Append the hPLZF transgene cDNA as one extra "transcript" record so reads
//      from the overexpression construct are captured by pseudoalignment instead
//      of being lost or misassigned to endogenous mouse Zbtb16 (docs/adr/0004).
//   3. Emit a tx2gene.tsv (transcript_id, gene_id, gene_name) covering every
//      mm10 transcript plus one synthetic hPLZF -> hPLZF row, for tximport's
//      transcript-to-gene summarization downstream.

process BUILD_TRANSCRIPTOME {
    tag 'mm10_vM25_plus_hPLZF'
    // gffread streams through the genome/GTF rather than loading an aligner-style
    // in-memory index (unlike STAR genome generation, which does need the
    // `big_mem` label in the sibling fork) -- default base.config resources
    // (8GB, retried up to 24GB) are enough headroom without over-requesting on
    // a modest machine.

    conda "${moduleDir}/../envs/gffread.yml"

    publishDir "${params.outdir}/salmon_kallisto/reference", mode: 'copy'

    input:
    path mm10_fasta
    path mm10_gtf
    path hplzf_fasta

    output:
    path 'transcriptome.fa', emit: transcriptome_fasta
    path 'tx2gene.tsv',      emit: tx2gene

    script:
    """
    # 1. Spliced transcript FASTA for all mm10 + GENCODE vM25 transcripts.
    gffread -w mm10_transcripts.fa -g ${mm10_fasta} ${mm10_gtf}

    # 2. Append the hPLZF transgene as one extra transcript entry. hPLZF.fa's
    #    header is already `>hPLZF` (single record) so no renaming is needed.
    cat mm10_transcripts.fa ${hplzf_fasta} > transcriptome.fa

    # 3. tx2gene map: one row per mm10 transcript (transcript_id, gene_id,
    #    gene_name pulled from GTF attributes of each `transcript` feature line),
    #    plus one synthetic row for hPLZF -> hPLZF (transcript IS the gene here,
    #    there is no separate mouse gene model for this construct).
    gawk 'BEGIN{FS="\\t"; OFS="\\t"}
        \$3 == "transcript" {
            attrs = \$9
            tid = ""; gid = ""; gname = ""
            if (match(attrs, /transcript_id "[^"]+"/)) {
                tid = substr(attrs, RSTART, RLENGTH)
                gsub(/transcript_id "|"/, "", tid)
            }
            if (match(attrs, /gene_id "[^"]+"/)) {
                gid = substr(attrs, RSTART, RLENGTH)
                gsub(/gene_id "|"/, "", gid)
            }
            if (match(attrs, /gene_name "[^"]+"/)) {
                gname = substr(attrs, RSTART, RLENGTH)
                gsub(/gene_name "|"/, "", gname)
            }
            if (gname == "") gname = gid
            if (tid != "") print tid, gid, gname
        }' ${mm10_gtf} > tx2gene_mm10.tsv

    printf 'hPLZF\\thPLZF\\thPLZF\\n' > tx2gene_hplzf.tsv
    cat tx2gene_mm10.tsv tx2gene_hplzf.tsv > tx2gene.tsv
    """

    stub:
    """
    printf '>ENSMUST00000stub.1\\nACGTACGTACGTACGTACGT\\n>hPLZF\\nACGTACGTACGTACGTACGT\\n' > transcriptome.fa
    printf 'ENSMUST00000stub.1\\tENSMUSG00000stub.1\\tStubGene\\nhPLZF\\thPLZF\\thPLZF\\n' > tx2gene.tsv
    """
}
