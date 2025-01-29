// 4. Transcriptome assembly with StringTie
process assembly_transcriptome_hisat2_stringtie {

  tag "Hisat2/StringTie - short reads"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work"
  publishDir "$params.outdir/evidence_data/transcriptomes/StringTie/short_reads/HISAT2"
  cpus 4

  input:
    val(bam_file)

  output:
    file("*_transcriptome.gtf"), emit: hisat2_stringtie_transcriptome_gtf
    file("*_transcriptome.AltCommands.gtf"), emit: hisat2_stringtie_alt_commands_gtf

  script:
    """
    /scripts/Stringtie.sh -t ${task.cpus} -o ${bam_file}_transcriptome.gtf -b $bam_file -r short
    /scripts/Stringtie_AltCommands.sh -t ${task.cpus} -o ${bam_file}_transcriptome.AltCommands.gtf -b $bam_file -r short
    """
}
