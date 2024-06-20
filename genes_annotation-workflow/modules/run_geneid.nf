// GeneID is a prediction tool that can't be trained easily so if the tool wasn't already ...
// ... trained for the organism the pipeline will annotate, this step will be skipped

process run_geneid {

  tag "Executing GeneID prediction"
  container 'avelt/geneid:latest'
  containerOptions "--volume ${genome_path}:/genome_path"
  publishDir "$params.outdir/geneid"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    val(geneid_param_file)

  output:
    file("geneid_predictions.gff3")

  script:
    """
    geneid /genome_path/$genome -P $geneid_param_file -3 > geneid_predictions.gff3
    """
}
