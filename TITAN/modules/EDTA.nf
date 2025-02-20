process EDTA {

  tag "Executing EDTA TE annotation on the following genome: $genome"
  container 'avelt/edta:latest'
  containerOptions "--memory=50g --volume ${genome_path}:/genome_path --volume ${projectDir}/scripts/:/scripts"
  cpus 5
  publishDir "${params.output_dir}/tmp", mode: 'copy'

  input:
    val(genome_path)
    val(genome)

  output:
    path("*TElib.fa"), emit : TElib_fasta
    path("*TEanno.gff3"), emit : TE_annotations_gff3
    path("*MAKER.masked"), emit : masked_genome

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running EDTA on $genome"
    CMD="/scripts/edta.sh -g /genome_path/$genome -n ${task.cpus}"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/edta.sh -g /genome_path/$genome -n ${task.cpus}
    cp *MAKER.masked ${params.output_dir}/assembly_masked.EDTA.fasta
    """
}
