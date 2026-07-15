// Builds the hybrid mm10 + GENCODE vM25 + hPLZF-transgene reference used by
// every other process in this fork. See docs/adr/0004.
//
// The hPLZF contig length is measured at run time (via `samtools faidx`)
// rather than hardcoded, so the synthetic GTF entry's end coordinate always
// matches whatever hPLZF.fa actually contains.

process BUILD_REFERENCE {
    tag "combined_genome"
    label 'alignment'
    label 'big_mem'

    conda "${moduleDir}/../envs/build_reference.yml"

    publishDir "${params.outdir}/reference", mode: 'copy', saveAs: { fn -> fn == 'star_index' ? null : fn }

    input:
    path mm10_fasta
    path mm10_gtf
    path hplzf_fasta

    output:
    path 'combined_genome.fa',       emit: fasta
    path 'combined_annotation.gtf',  emit: gtf
    path 'star_index',               emit: index

    script:
    """
    set -euo pipefail

    # 1) Combined genome FASTA: mm10 main contigs (chr1-19,X,Y,M) + the hPLZF
    #    transgene as its own contig.
    cat ${mm10_fasta} ${hplzf_fasta} > combined_genome.fa

    # 2) Measure the hPLZF contig length/name so the synthetic GTF entry is
    #    generated, not hand-typed.
    samtools faidx ${hplzf_fasta}
    hplzf_id=\$(cut -f1 ${hplzf_fasta}.fai | head -n1)
    hplzf_len=\$(cut -f2 ${hplzf_fasta}.fai | head -n1)

    # 3) Combined GTF: GENCODE vM25 annotation + one synthetic gene/transcript/
    #    exon triplet for hPLZF, spanning the whole contig as a single exon
    #    (gene_biotype "transgene", matching the Plasmidsaurus matrix's
    #    hPLZF row semantics -- see docs/adr/0004).
    cp ${mm10_gtf} combined_annotation.gtf
    {
      printf '%s\\tcustom\\tgene\\t1\\t%s\\t.\\t+\\t.\\tgene_id "hPLZF"; gene_name "hPLZF"; gene_biotype "transgene";\\n' "\${hplzf_id}" "\${hplzf_len}"
      printf '%s\\tcustom\\ttranscript\\t1\\t%s\\t.\\t+\\t.\\tgene_id "hPLZF"; transcript_id "hPLZF.1"; gene_name "hPLZF"; gene_biotype "transgene"; transcript_biotype "transgene";\\n' "\${hplzf_id}" "\${hplzf_len}"
      printf '%s\\tcustom\\texon\\t1\\t%s\\t.\\t+\\t.\\tgene_id "hPLZF"; transcript_id "hPLZF.1"; gene_name "hPLZF"; gene_biotype "transgene"; transcript_biotype "transgene"; exon_number 1; exon_id "hPLZF.1.1";\\n' "\${hplzf_id}" "\${hplzf_len}"
    } >> combined_annotation.gtf

    # 4) STAR genome index over the combined reference.
    mkdir -p star_index
    STAR \\
      --runMode genomeGenerate \\
      --genomeDir star_index \\
      --genomeFastaFiles combined_genome.fa \\
      --sjdbGTFfile combined_annotation.gtf \\
      --sjdbOverhang ${params.star_sjdb_overhang} \\
      --runThreadN ${task.cpus}
    """

    stub:
    """
    mkdir -p star_index
    printf '>hPLZF\\nACGTACGTACGT\\n' > combined_genome.fa
    printf 'hPLZF\\tcustom\\tgene\\t1\\t12\\t.\\t+\\t.\\tgene_id "hPLZF"; gene_name "hPLZF"; gene_biotype "transgene";\\n' > combined_annotation.gtf
    touch star_index/SA star_index/SAindex star_index/Genome star_index/genomeParameters.txt
    """
}
