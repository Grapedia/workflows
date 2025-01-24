process gffread_convert_gff3_to_cds_fasta {

  tag "Executing gffread to convert liftoff gff3 to cds.fasta"
  container 'quay.io/biocontainers/gffread:0.12.7--h077b44d_6'
  containerOptions "--volume $genome_path:/genome_path --volume ${projectDir}/scripts:/scripts"
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
    liftoff_gff3_path=\$(/scripts/retrieve_path_liftoff_gff3.sh ${liftoff_gff3})
    gffread -x ${genome}.CDS.fasta -g /genome_path/${genome} \${liftoff_gff3_path}
    gzip ${genome}.CDS.fasta
    """
}
