// 4. Transcriptome assembly with StringTie
process assembly_transcriptome_minimap2_stringtie {

  tag "Minimap2/StringTie - long reads"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work"
  publishDir "$params.outdir/evidence_data/transcriptomes/StringTie/long_reads"
  cpus 4

  input:
    val(bam_file)

  output:
    path("*_transcriptome.gtf", emit: minimap2_stringtie_transcriptome_gtf)
    path("*_transcriptome.AltCommands.gtf", emit: minimap2_stringtie_alt_commands_gtf)

  when:
  has_long_reads

  script:
    """
    /scripts/Stringtie.sh -t ${task.cpus} -o ${bam_file}_transcriptome.gtf -b $bam_file -r long
    /scripts/Stringtie_AltCommands.sh -o ${bam_file}_transcriptome.AltCommands.gtf -b $bam_file -r long
    """
}
