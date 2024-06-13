// 4. Transcriptome assembly
process assembly_transcriptome_unstranded {

  tag "psiclass transcriptome assembly - unstranded"
  container 'quay.io/biocontainers/psiclass:1.0.2--h87f3376_2'
  containerOptions "--volume $params.outdir/evidence_data/RNAseq_unstranded/alignments/new_assembly:/alignments"
  publishDir "$params.outdir/evidence_data/transcriptomes/RNAseq_unstranded"
  cpus 4

  input:
    tuple val(sample_ID), val(stranded_or_unstranded), val(paired_or_single), path(bam)

  output:
    file("RNAseq_unstranded_vote.gtf")

  script:
    """
    bam="\$(ls -1 /alignments/*bam | tr '\n' ',' | sed 's/,\$//')"
    psiclass -p ${task.cpus} -b \${bam} -o RNAseq_unstranded
    """
}
