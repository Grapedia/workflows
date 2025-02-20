// 4. Transcriptome assembly with PsiCLASS
process assembly_transcriptome_star_psiclass {

  tag "STAR/PsiCLASS - short reads - on ${sample_ID}"
  container 'avelt/psiclass_samtools:latest'
  publishDir { 
    if (strand_type == "unstranded") {
      return "$projectDir/FINAL_OUTPUT/transcriptomes/STAR_PsiCLASS/unstranded"
    } else if (strand_type in ["stranded_forward", "stranded_reverse"]) {
      return "$projectDir/FINAL_OUTPUT/transcriptomes/STAR_PsiCLASS/stranded"
    }
  }
  cpus 4

  input:
    tuple val(sample_ID), path(bam_file), val(strand_type)

  output:
    tuple val(sample_ID), file("${sample_ID}_vote.gtf"), val(strand_type), emit: psiclass_assemblies

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running PsiCLASS transcriptome assembly on $sample_ID"
    CMD="/PsiCLASS-1.0.2/psiclass -p ${task.cpus} -b ${bam_file} -o ${sample_ID}"
    echo "[\$DATE] Executing: \$CMD"
    /PsiCLASS-1.0.2/psiclass -p ${task.cpus} -b ${bam_file} -o ${sample_ID}
    """
}
