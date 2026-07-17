process validate_final_annotation {
  label 'process_low'
  tag "TITAN final annotation validation"

  container params.container_python
  publishDir "${params.output_dir}/validation", mode: 'copy'

  input:
    path(genome)
    path(annotation)
    path(proteins_all)
    path(proteins_main)
    path(validation_script)

  output:
    path "final_annotation_validation.json", emit: json_report
    path "final_annotation_validation.txt", emit: text_report
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    python3 ${validation_script} \\
      --genome ${genome} \\
      --annotation ${annotation} \\
      --proteins-all ${proteins_all} \\
      --proteins-main ${proteins_main} \\
      --json-report final_annotation_validation.json \\
      --text-report final_annotation_validation.txt
    script_sha256=\$(sha256sum ${validation_script} | awk '{print \$1}')
    printf '"%s":\n  container: "%s"\n  validation_script: "%s"\n  validation_script_sha256: "%s"\n' \\
      "${task.process}" "${task.container}" "${validation_script}" "\${script_sha256}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    python3 ${validation_script} \\
      --genome ${genome} \\
      --annotation ${annotation} \\
      --proteins-all ${proteins_all} \\
      --proteins-main ${proteins_main} \\
      --json-report final_annotation_validation.json \\
      --text-report final_annotation_validation.txt
    script_sha256=\$(sha256sum ${validation_script} | awk '{print \$1}')
    printf '"%s":\n  container: "%s"\n  validation_script: "%s"\n  validation_script_sha256: "%s"\n' \\
      "${task.process}" "${task.container}" "${validation_script}" "\${script_sha256}" > versions.yml
    """
}
