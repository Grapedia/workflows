process glimmerhmm_training {

  tag "Executing GlimmerHMM training"
  container 'avelt/glimmerhmm_gffutils:latest'
  containerOptions "--volume ${projectDir}/work:/work --volume $params.outdir:/outdir --volume $genome_path:/genome_path --volume ${projectDir}/scripts/:/scripts --volume $params.outdir/evidence_data/transcriptomes/rnaseq_stranded:/transcriptomes_rnaseq_stranded"
  publishDir "$params.outdir/GlimmerHMM/"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    val(stranded_gff3)

  output:
    path("training")

  script:
    """
    /scripts/trainglimmerhmm.sh -a /genome_path/$genome -t /transcriptomes_rnaseq_stranded/transcriptome_RNAseq_stranded.gff3 -d training -i /outdir
    """
}
