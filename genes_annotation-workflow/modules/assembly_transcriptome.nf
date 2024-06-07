// 4. Transcriptome assembly
process assembly_transcriptome {

  tag "psiclass transcriptome assembly"
  container 'quay.io/biocontainers/psiclass:1.0.2--h87f3376_2'
  containerOptions "--volume $params.outdir/evidence_data/RNAseq_$stranded_or_unstranded/alignments/new_assembly:/alignments"
  publishDir "$params.outdir/evidence_data/transcriptomes/$stranded_or_unstranded"
  cpus 4

  input:

  output:
    file("$stranded_or_unstranded_vote.gtf")

  script:
    """
    psiclass -p ${task.cpus} -b -o {params.prefix_path}
    """
}
