process agat_convert_gff3_to_cds_fasta {
  label 'process_low'

  tag "Executing agat to convert liftoff gff3 to cds.fasta"
  container params.container_agat
  publishDir "${params.output_dir}/intermediate_files/liftoff/gff3_to_cds_fasta", mode: "copy", enabled: params.publish_intermediates
  input:
    path(genome)
    path(liftoff_gff3)

  output:
    path("${genome.simpleName}.CDS.fasta.gz"), emit: cds_fasta
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running agat_convert_gff3_to_cds_fasta on ${genome}"

    fold -w 80 ${genome} > reformatted.fa
    CMD="agat_sp_extract_sequences.pl -g ${liftoff_gff3} -f reformatted.fa -o ${genome.simpleName}.CDS.fasta"
    echo "[\$DATE] Executing: \$CMD"
    agat_sp_extract_sequences.pl -g ${liftoff_gff3} -f reformatted.fa -o ${genome.simpleName}.CDS.fasta
    gzip ${genome.simpleName}.CDS.fasta
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf ">stub_cds\\nATGGCC\\n" | gzip -c > ${genome.simpleName}.CDS.fasta.gz
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
