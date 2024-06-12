process glimmerhmm_training {

  tag "Executing GlimmerHMM training"
  container 'avelt/glimmerhmm_gffutils:latest'
  containerOptions "--volume $genome_path:/genome_path"
  publishDir "$params.outdir/evidence_data/databases/"
  cpus 4

  input:
    val(genome_path)
    val(genome)

  output:
    path("new_assembly_database")

  script:
    """
    gmap_build --nthreads ${task.cpus} --dir \$PWD --genomedb new_assembly_database /genome_path/$genome
    """
}
