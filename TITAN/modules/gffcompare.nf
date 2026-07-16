// 2. gffcompare
process gffcompare {

  tag "gffcompare on STAR/PsiCLASS GTF to merge them"
  container params.container_gffcompare
  cpus 4

  publishDir "${params.output_dir}/tmp", mode: 'copy'
  publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename ->
    if (filename == 'stranded_merged_output.combined.gtf') {
      return 'merged_star_psiclass_stranded.gtf'
    }
    if (filename == 'unstranded_merged_output.combined.gtf') {
      return 'merged_star_psiclass_unstranded.gtf'
    }
    return null
  }

  input:
    path(stranded_gtfs)
    path(unstranded_gtfs)

  output:
    path "stranded_merged_output.combined.gtf", emit: star_psiclass_stranded
    path "unstranded_merged_output.combined.gtf", optional: true, emit: star_psiclass_unstranded

  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running gffcompare on STAR/PsiCLASS GTF to merge them"
    /gffcompare-0.12.6/gffcompare -o stranded_merged_output ${stranded_gtfs}

    if [[ -n "${unstranded_gtfs}" ]]; then
      /gffcompare-0.12.6/gffcompare -o unstranded_merged_output ${unstranded_gtfs}
    else
      : > unstranded_merged_output.combined.gtf
    fi
    """

  stub:
    """
    printf "chr1\\tgffcompare\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"psiclass_stranded_gene\\"; transcript_id \\"psiclass_stranded_tx\\";\\n" > stranded_merged_output.combined.gtf
    printf "chr1\\tgffcompare\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"psiclass_unstranded_gene\\"; transcript_id \\"psiclass_unstranded_tx\\";\\n" > unstranded_merged_output.combined.gtf
    """
}
