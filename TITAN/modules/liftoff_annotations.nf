// Previous annotations will be lifted over (mapped) on the new assembly. It takes in ...
// ... input:
//   * the reference assembly
//   * the target assembly
// * the annotations to lift over the target assembly
process liftoff_annotations {
  label 'process_low'

  cache 'deep'
  tag "Executing liftoff on the new assembly $genome"
  container params.container_liftoff
  publishDir "${params.output_dir}", mode: 'copy'

  input:
    path(genome_fasta)
    val(genome)
    path(previous_assembly_fasta)
    val(previous_assembly)
    path(previous_annotations_gff3)
    val(previous_annotations)

  output:
    path("liftoff_previous_annotations.gff3"), emit : liftoff_previous_annotations
    path("unmapped_features.txt"), emit : unmapped_features
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running liftoff on $genome"

    CMD="liftoff -g ${previous_annotations_gff3} -o liftoff_previous_annotations.gff3 -u unmapped_features.txt ${genome_fasta} ${previous_assembly_fasta}"

    echo "[\$DATE] Executing: \$CMD"
    liftoff -g ${previous_annotations_gff3} -o liftoff_previous_annotations.gff3 -u unmapped_features.txt ${genome_fasta} ${previous_assembly_fasta}
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml

    """

  stub:
    """
    printf "##gff-version 3\\nchr1\\tliftoff\\tgene\\t1\\t10\\t.\\t+\\t.\\tID=liftoff_stub_gene\\n" > liftoff_previous_annotations.gff3
    printf "stub_unmapped\\n" > unmapped_features.txt
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
