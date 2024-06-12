process glimmerhmm_training {

  tag "Executing GlimmerHMM training"
  container 'avelt/glimmerhmm_gffutils:latest'
  containerOptions "--volume $params.outdir:/outdir --volume $genome_path:/genome_path --volume ${projectDir}/scripts/:/scripts -volume $params.outdir/evidence_data/transcriptomes/RNAseq_stranded:/transcriptomes_RNAseq_stranded"
  publishDir "$params.outdir/GlimmerHMM/"
  cpus 4

  input:
    val(genome_path)
    val(genome)

  output:
    path("training")

  script:
    """
    /scripts/trainglimmerhmm.sh -a /genome_path/$genome -t /transcriptomes_RNAseq_stranded/*gff3 -d training -i /outdir
    """
}
