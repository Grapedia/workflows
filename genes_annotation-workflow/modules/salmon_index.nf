process salmon_index {

  tag "Executing salmon indexing on $cds_fasta"
  container 'quay.io/biocontainers/salmon:1.10.3--haf24da9_3'
  publishDir "$params.outdir/salmon_index/"
  cpus 4

  input:
    val(cds_fasta)

  output:
    path("salmon_index")

  script:
    """
    salmon index -t $cds_fasta -i salmon_index
    """
}
