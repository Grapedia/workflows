// 1. Building database for reference genome with gmap_build
process gmap_build_database {

  tag "gmap_build on $genome"
  container 'quay.io/biocontainers/gmap:2020.10.14--pl526h2f06484_0'
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
