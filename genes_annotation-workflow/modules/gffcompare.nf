// 2. gffcompare
process gffcompare {

  tag "gffcompare on STAR/PsiCLASS GTF to merge them"
  container 'avelt/gffcompare:latest'
  containerOptions "--volume ${projectDir}/work:/work --volume ${projectDir}/scripts:/scripts --volume $projectDir/FINAL_OUTPUT/transcriptomes/STAR_PsiCLASS/stranded/:/STAR_PsiCLASS_stranded --volume $projectDir/FINAL_OUTPUT/transcriptomes/STAR_PsiCLASS/unstranded/:/STAR_PsiCLASS_unstranded"
  publishDir "$projectDir/FINAL_OUTPUT/transcriptomes/STAR_PsiCLASS/"
  cpus 4

  input:
    val(merged_gtf)

  output:
    file("*merged_output*")

  script:
    """
    gtf_files=\$(/scripts/retrieve_path_gffcompare.sh /STAR_PsiCLASS_stranded/)
    /gffcompare-0.12.6/gffcompare -o stranded_merged_output \${gtf_files}
    gtf_files=\$(/scripts/retrieve_path_gffcompare.sh /STAR_PsiCLASS_unstranded/)
    /gffcompare-0.12.6/gffcompare -o unstranded_merged_output \${gtf_files}
    """
}
