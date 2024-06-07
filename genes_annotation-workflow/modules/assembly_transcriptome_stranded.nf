// 4. Transcriptome assembly
process assembly_transcriptome_stranded {

  tag "psiclass transcriptome assembly"
  container 'quay.io/biocontainers/psiclass:1.0.2--h87f3376_2'
  containerOptions "--volume $params.outdir/evidence_data/RNAseq_stranded/alignments/new_assembly:/alignments"
  publishDir "$params.outdir/evidence_data/transcriptomes/RNAseq_stranded"
  cpus 4

  input:
    tuple val(sample_ID), val(stranded_or_unstranded), val(paired_or_single), path(bam)

  output:
    file("RNAseq_stranded_vote.gtf")

  script:
    """
    psiclass -p ${task.cpus} -b `ls -1 /alignments/*bam | tr '\n' ',' | sed 's/,\$//'` -o RNAseq_stranded
    """
}
