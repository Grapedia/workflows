// 4. Transcriptome assembly with PsiCLASS
process assembly_transcriptome_star_psiclass {
  label 'process_transcriptome'

  tag "STAR/PsiCLASS - short reads - on ${sample_ID}"
  container params.container_psiclass_samtools
  publishDir {
    if (strand_type == "unstranded") {
      return "${params.output_dir}/intermediate_files/evidence_data/transcriptomes/STAR_PsiCLASS/unstranded"
    } else if (strand_type in ["stranded_forward", "stranded_reverse"]) {
      return "${params.output_dir}/intermediate_files/evidence_data/transcriptomes/STAR_PsiCLASS/stranded"
    }
  }, mode: "copy", enabled: params.publish_intermediates
  input:
    tuple val(sample_ID), path(bam_file), val(strand_type)
    val(psiclass_vd_option)
    val(psiclass_c_option)

  output:
    tuple val(sample_ID), path("${sample_ID}_vote.gtf"), val(strand_type), emit: psiclass_assemblies
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running PsiCLASS transcriptome assembly on $sample_ID"
    CMD="/PsiCLASS-1.0.2/psiclass -p ${task.cpus} -b ${bam_file} -o ${sample_ID} --vd ${psiclass_vd_option} -c ${psiclass_c_option} --primaryParalog"
    echo "[\$DATE] Executing: \$CMD"
    psiclass_status=0
    /PsiCLASS-1.0.2/psiclass -p ${task.cpus} -b ${bam_file} -o ${sample_ID} --vd ${psiclass_vd_option} -c ${psiclass_c_option} --primaryParalog || psiclass_status=\$?
    if [[ "\${psiclass_status}" -ne 0 && ! -s "${sample_ID}_vote.gtf" ]]; then
      trusted_splice_count=\$(find splice -type f -name '*trusted_splice' -size +0c 2>/dev/null | wc -l)
      if [[ "\${trusted_splice_count}" -eq 0 ]]; then
        echo "[\$DATE] WARNING: PsiCLASS exited \${psiclass_status} but no trusted splice junctions were produced; emitting an empty GTF for $sample_ID." >&2
        : > ${sample_ID}_vote.gtf
      else
        exit "\${psiclass_status}"
      fi
    fi
    printf '"%s":\n  psiclass: "1.0.2"\n  psiclass_vd_option: "%s"\n  psiclass_c_option: "%s"\n  psiclass_exit_status: "%s"\n' "${task.process}" "${psiclass_vd_option}" "${psiclass_c_option}" "\${psiclass_status}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "chr1\\tPsiCLASS\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_psiclass_gene\\"; transcript_id \\"${sample_ID}_psiclass_tx\\";\\n" > ${sample_ID}_vote.gtf
    printf '"%s":\n  psiclass: "stub"\n  psiclass_vd_option: "%s"\n  psiclass_c_option: "%s"\n' "${task.process}" "${psiclass_vd_option}" "${psiclass_c_option}" > versions.yml
    """
}
