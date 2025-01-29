// 4. Transcriptome assembly with StringTie
process assembly_transcriptome_star_stringtie {

  tag "STAR/StringTie - short reads"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work"
  publishDir "$params.outdir/evidence_data/transcriptomes/StringTie/short_reads/STAR"
  cpus 4

  input:
    val(bam_file)

  output:
    file("*_transcriptome.gtf"), emit: star_stringtie_transcriptome_gtf
    file("*_transcriptome.AltCommands.gtf"), emit: star_stringtie_alt_commands_gtf

  script:
    """
    /scripts/Stringtie.sh -t ${task.cpus} -o ${bam_file}_transcriptome.gtf -b $bam_file -r short
    /scripts/Stringtie_AltCommands.sh -t ${task.cpus} -o ${bam_file}_transcriptome.AltCommands.gtf -b $bam_file -r short
    """
}
