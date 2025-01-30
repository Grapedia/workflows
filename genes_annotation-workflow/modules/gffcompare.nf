// 2. gffcompare
process gffcompare {

  tag "gffcompare on STAR/PsiCLASS GTF to merge them"
  container 'avelt/gffcompare:latest'
  containerOptions "--volume ${projectDir}/work:/work --volume ${projectDir}/scripts:/scripts --volume $projectDir/FINAL_OUTPUT/transcriptomes/STAR_PsiCLASS/stranded/:/STAR_PsiCLASS/stranded/ --volume $projectDir/FINAL_OUTPUT/transcriptomes/STAR_PsiCLASS/unstranded/:/STAR_PsiCLASS/unstranded/"
  publishDir "$projectDir/FINAL_OUTPUT/transcriptomes/STAR_PsiCLASS/"
  cpus 4

  input:
    val(merged_gtf)

  output:
    file("*merged_output*")

  script:
    """
    /gffcompare-0.12.6/gffcompare -o stranded_merged_output /STAR_PsiCLASS/stranded/*gtf
    /gffcompare-0.12.6/gffcompare -o unstranded_merged_output /STAR_PsiCLASS/unstranded/*gtf
    """
}
