// 4. Transcriptome assembly with StringTie
process assembly_transcriptome_star_stringtie {

  tag "StringTie transcriptome assembly with STAR alignements"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work"
  publishDir "$params.outdir/evidence_data/transcriptomes/StringTie"
  cpus 4

  input:
    val(bam_file)

  output:
    file("*_transcriptome.gtf")

  script:
    """
    /scripts/Stringtie.sh -t ${task.cpus} -o ${bam_file}_transcriptome.gtf -b $bam_file
    """
}
