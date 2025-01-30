// Previous annotations will be lifted over (mapped) on the new assembly. It takes in ...
// ... input:
//   * the reference assembly
//   * the target assembly
// * the annotations to lift over the target assembly
process liftoff_annotations {

  cache 'deep'
  tag "Executing liftoff on the new assembly $genome"
  container 'quay.io/biocontainers/liftoff:1.5.1--py_0'
  containerOptions "--volume $genome_path:/genome_path --volume $previous_assembly_path:/previous_assembly_path --volume $previous_annotations_path:/previous_annotations_path"
  publishDir "$projectDir/FINAL_OUTPUT/liftoff/"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    val(previous_assembly_path)
    val(previous_assembly)
    val(previous_annotations_path)
    val(previous_annotations)

  output:
    path("liftoff_previous_annotations.gff3"), emit : liftoff_previous_annotations
    path("unmapped_features.txt"), emit : unmapped_features

  script:
    """
    liftoff -g /previous_annotations_path/$previous_annotations -o liftoff_previous_annotations.gff3 -u unmapped_features.txt /genome_path/$genome /previous_assembly_path/$previous_assembly
    """
}
