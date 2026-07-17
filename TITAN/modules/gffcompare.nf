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
    path(stranded_gtfs, stageAs: "stranded_gtfs/*")
    path(unstranded_gtfs, stageAs: "unstranded_gtfs/*")

  output:
    path "stranded_merged_output.combined.gtf", emit: star_psiclass_stranded
    path "unstranded_merged_output.combined.gtf", emit: star_psiclass_unstranded
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running gffcompare on STAR/PsiCLASS GTF to merge them"

    write_nonempty_gtf_list() {
      local input_dir="\$1"
      local output_list="\$2"
      local required="\$3"
      local count=0

      : > "\${output_list}"
      while IFS= read -r -d '' gtf_file; do
        if [[ -s "\${gtf_file}" ]]; then
          printf '%s\\n' "\${gtf_file}" >> "\${output_list}"
          count=\$((count + 1))
        fi
      done < <(find "\${input_dir}" -type f -name '*.gtf' -print0 | sort -z)

      if [[ "\${required}" == "true" && "\${count}" -eq 0 ]]; then
        echo "[\$DATE] ERROR: no non-empty GTF files found in \${input_dir}" >&2
        exit 1
      fi
      [[ "\${count}" -gt 0 ]]
    }

    write_nonempty_gtf_list stranded_gtfs stranded_gtfs.txt true
    /gffcompare-0.12.6/gffcompare -o stranded_merged_output -i stranded_gtfs.txt

    if write_nonempty_gtf_list unstranded_gtfs unstranded_gtfs.txt false; then
      /gffcompare-0.12.6/gffcompare -o unstranded_merged_output -i unstranded_gtfs.txt
    else
      : > unstranded_merged_output.combined.gtf
    fi
    printf '"%s":\n  gffcompare: "0.12.6"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    : > stranded_merged_output.combined.gtf
    : > unstranded_merged_output.combined.gtf
    printf '"%s":\n  gffcompare: "stub"\n' "${task.process}" > versions.yml
    """
}
