process RepeatMasker {

  tag "Executing RepeatMasker on $genome"
  container 'quay.io/biocontainers/repeatmasker:4.1.0--pl526_0'
  containerOptions "--volume ${genome_path}:/genome_path"
  publishDir "$params.outdir/RepeatMasker/"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    val(TEs_lib_EDTA)

  output:
    file("*.masked"), emit : masked_genome
    file("*.gff"), emit : RepeatMasker_annot

  script:
    """
    RepeatMasker -pa ${task.cpus} -gff -lib $TEs_lib_EDTA /genome_path/$genome
    """
}
