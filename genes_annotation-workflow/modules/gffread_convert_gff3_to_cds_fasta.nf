process gffread_convert_gff3_to_cds_fasta {

  tag "Executing gffread to convert liftoff gff3 to cds.fasta"
  container 'quay.io/biocontainers/gffread:0.12.7--h077b44d_6'
  containerOptions "--volume $genome_path:/genome_path --volume ${projectDir}/scripts:/scripts --volume ${projectDir}/work:/work"
  publishDir "$params.outdir/liftoff/gff3_to_cds_fasta"
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
    echo "[\$DATE] Running gffread_convert_gff3_to_cds_fasta on $genome" >> ${params.logfile} 2>&1

    liftoff_gff3_path=\$(/scripts/retrieve_path.sh ${liftoff_gff3})
    CMD="gffread -x ${genome}.CDS.fasta -g /genome_path/${genome} \${liftoff_gff3_path}"
    echo "[\$DATE] Executing: \$CMD" >> ${params.logfile} 2>&1
    gffread -x ${genome}.CDS.fasta -g /genome_path/${genome} \${liftoff_gff3_path}
    gzip ${genome}.CDS.fasta
    """
}
