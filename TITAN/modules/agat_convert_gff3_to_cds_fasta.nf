process agat_convert_gff3_to_cds_fasta {

  tag "Executing agat to convert liftoff gff3 to cds.fasta"
  container 'quay.io/biocontainers/agat:1.2.0--pl5321hdfd78af_0'
  containerOptions "--volume $genome_path:/genome_path --volume ${projectDir}/scripts:/scripts --volume ${projectDir}/work:/work"
  publishDir "${params.output_dir}/intermediate_files/liftoff/gff3_to_cds_fasta"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    val(liftoff_gff3)

  output:
    file("${genome}.CDS.fasta.gz")

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running agat_convert_gff3_to_cds_fasta on $genome"

    liftoff_gff3_path=\$(/scripts/retrieve_path.sh ${liftoff_gff3})
    fold -w 80 /genome_path/${genome} > reformatted.fa
    CMD="agat_sp_extract_sequences.pl -g \${liftoff_gff3_path} -f reformatted.fa -o ${genome}.CDS.fasta"
    echo "[\$DATE] Executing: \$CMD"
    agat_sp_extract_sequences.pl -g \${liftoff_gff3_path} -f reformatted.fa -o ${genome}.CDS.fasta
    gzip ${genome}.CDS.fasta
    """
}