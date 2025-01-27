process salmon_index {

  tag "Executing salmon indexing on $cds_fasta"
  container 'quay.io/biocontainers/salmon:1.10.3--haf24da9_3'
  containerOptions "--volume ${projectDir}/work:/work --volume ${projectDir}/scripts:/scripts"
  publishDir "$params.outdir/salmon_index/"
  cpus 4

  input:
    val(cds_fasta)

  output:
    path("salmon_index")

  script:
    """
    index_path=\$(/scripts/retrieve_path.sh ${cds_fasta})
    salmon index -t \${index_path} -i salmon_index
    """
}
