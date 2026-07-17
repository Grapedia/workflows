process agat_convert_gff3_to_cds_fasta {
  label 'process_low'

  tag "Executing agat to convert liftoff gff3 to cds.fasta"
  container params.container_agat
  publishDir "${params.output_dir}/intermediate_files/liftoff/gff3_to_cds_fasta", mode: "copy", enabled: params.publish_intermediates
  input:
    path(genome)
    path(liftoff_gff3)
    path(clean_liftoff_gff3_script)

  output:
    path("${genome.simpleName}.CDS.fasta.gz"), emit: cds_fasta
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running agat_convert_gff3_to_cds_fasta on ${genome}"

    command -v python3 >/dev/null || { echo "python3 is required in the AGAT container to clean Liftoff GFF3 before AGAT extraction" >&2; exit 127; }
    python3 ${clean_liftoff_gff3_script} \\
      --input ${liftoff_gff3} \\
      --output cleaned.OK.gff3 \\
      --removed-ids removed_feature_ids.txt

    fold -w 80 ${genome} > reformatted.fa
    CMD="agat_sp_extract_sequences.pl -g cleaned.OK.gff3 -f reformatted.fa -o ${genome.simpleName}.CDS.fasta"
    echo "[\$DATE] Executing: \$CMD"
    agat_sp_extract_sequences.pl -g cleaned.OK.gff3 -f reformatted.fa -o ${genome.simpleName}.CDS.fasta
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
