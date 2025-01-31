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
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running gffcompare on STAR/PsiCLASS GTF to merge them" >> ${params.logfile} 2>&1
    gtf_files=\$(/scripts/retrieve_path_gffcompare.sh /STAR_PsiCLASS_stranded/)
    CMD="/gffcompare-0.12.6/gffcompare -o stranded_merged_output \${gtf_files}"
    echo "[\$DATE] Executing: \$CMD" >> ${params.logfile} 2>&1
    /gffcompare-0.12.6/gffcompare -o stranded_merged_output \${gtf_files}

    if [ -d "/STAR_PsiCLASS_unstranded/" ] && [ "\$(ls -A /STAR_PsiCLASS_unstranded/ 2>/dev/null)" ]; then    
      gtf_files=\$(/scripts/retrieve_path_gffcompare.sh /STAR_PsiCLASS_unstranded/)
      CMD="/gffcompare-0.12.6/gffcompare -o unstranded_merged_output \${gtf_files}"
      echo "[\$DATE] Executing: \$CMD" >> ${params.logfile} 2>&1
      /gffcompare-0.12.6/gffcompare -o unstranded_merged_output \${gtf_files}
    fi
    """
}
