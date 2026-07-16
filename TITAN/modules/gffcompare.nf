// 2. gffcompare
process gffcompare {
  label 'process_merge'

  tag "gffcompare STAR/PsiCLASS: stranded and optional unstranded GTFs"
  container params.container_gffcompare
  publishDir "${params.output_dir}/tmp", mode: 'copy', enabled: params.publish_intermediates
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
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running gffcompare on STAR/PsiCLASS GTF to merge them"
    /gffcompare-0.12.6/gffcompare -o stranded_merged_output ${stranded_gtfs}

    unstranded_files=( ${unstranded_gtfs} )
    if [[ \${#unstranded_files[@]} -gt 0 && -s "\${unstranded_files[0]}" ]]; then
      /gffcompare-0.12.6/gffcompare -o unstranded_merged_output ${unstranded_gtfs}
    else
      : > unstranded_merged_output.combined.gtf
    fi
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    printf "chr1\\tgffcompare\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"psiclass_stranded_gene\\"; transcript_id \\"psiclass_stranded_tx\\";\\n" > stranded_merged_output.combined.gtf
    printf "chr1\\tgffcompare\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"psiclass_unstranded_gene\\"; transcript_id \\"psiclass_unstranded_tx\\";\\n" > unstranded_merged_output.combined.gtf
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
