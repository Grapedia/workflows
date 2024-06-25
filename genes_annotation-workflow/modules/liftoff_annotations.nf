// Previous annotations will be lifted over (mapped) on the new assembly. It takes in ...
// ... input:
//   * the reference assembly
//   * the target assembly
// * the annotations to lift over the target assembly
process liftoff_annotations {

  tag "Executing liftoff on the new assembly"
  container 'quay.io/biocontainers/liftoff:1.5.1--py_0'
  containerOptions "--volume $genome_path:/genome_path --volume $annotations_path:/annotations_path"
  publishDir "$params.outdir/liftoff/"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    val(annotations_path)
    val(previous_assembly)
    val(previous_annotations)

  output:
    file("liftoff_previous_annotations.gff3"), emit : liftoff_previous_annotations
    file("unmapped_features.txt"), emit : unmapped_features

  script:
    """
    liftoff -g /annotations_path/$previous_annotations -o liftoff_previous_annotations.gff3 -u unmapped_features.txt /genome_path/$genome /genome_path/$previous_assembly
    """
}
