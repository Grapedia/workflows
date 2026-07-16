process validate_final_annotation {
  label 'process_low'
  tag "TITAN final annotation validation"

  publishDir "${params.output_dir}/validation", mode: 'copy'

  input:
    path(genome)
    path(annotation)
    path(proteins_all)
    path(proteins_main)

  output:
    path "final_annotation_validation.json", emit: json_report
    path "final_annotation_validation.txt", emit: text_report
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    python3 ${projectDir}/scripts/validate_final_annotation.py \\
      --genome ${genome} \\
      --annotation ${annotation} \\
      --proteins-all ${proteins_all} \\
      --proteins-main ${proteins_main} \\
      --json-report final_annotation_validation.json \\
      --text-report final_annotation_validation.txt
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    python3 ${projectDir}/scripts/validate_final_annotation.py \\
      --genome ${genome} \\
      --annotation ${annotation} \\
      --proteins-all ${proteins_all} \\
      --proteins-main ${proteins_main} \\
      --json-report final_annotation_validation.json \\
      --text-report final_annotation_validation.txt
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
