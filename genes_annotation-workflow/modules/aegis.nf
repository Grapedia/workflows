// Final step : gff3 merging with Aegis to create final GFF3 file
process aegis {

  tag "Generation of final GFF3 file with Aegis"
  // this image contains Aegis, Diamond and gffread
  container 'avelt/aegis:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work"
  publishDir "$projectDir/FINAL_OUTPUT/final_GFF3_file/"
  cpus 4

  input:
    tuple val(sample_ID), path(bam_file), val(strand_type)

  output:
    tuple val(sample_ID), path("*_transcriptome.gtf"), val(strand_type), emit: hisat2_stringtie_transcriptome_gtf
    tuple val(sample_ID), path("*_transcriptome.AltCommands.gtf"), val(strand_type), emit: hisat2_stringtie_alt_commands_gtf

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running Aegis to create final GFF3 file" >> ${params.logfile} 2>&1
    # convert all GTF file in GFF3 format with gffread -E {gtf_file} -o- > {gff3_file}
    /scripts/Aegis1.py -h
    /scripts/Aegis2.sh -h
    /scripts/Aegis3.py -h
    """
}
