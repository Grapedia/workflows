// Transcriptomes merging with StringTie
process Stringtie_merging_long_reads {
  label 'process_merge'

  tag "Minimap2/StringTie merge: long-read default and alt GTFs"
  container params.container_stringtie
  stageInMode 'copy'
  publishDir "${params.output_dir}/tmp", mode: 'copy', enabled: params.publish_intermediates
  publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename ->
    if (filename == 'merged_transcriptomes.minimap2.long_reads.default_args.gtf') {
      return 'merged_minimap2_stringtie_long_reads_default.gtf'
    }
    if (filename == 'merged_transcriptomes.minimap2.long_reads.alt_args.gtf') {
      return 'merged_minimap2_stringtie_long_reads_alt.gtf'
    }
    return null
  }

  input:
    path(default_gtfs, stageAs: "default_gtfs/*")
    path(alt_gtfs, stageAs: "alt_gtfs/*")

  output:
    path "merged_transcriptomes.minimap2.long_reads.default_args.gtf", emit: default_args_gtf
    path "merged_transcriptomes.minimap2.long_reads.alt_args.gtf", emit: alt_args_gtf
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running Minimap2/StringTie merging - long reads transcriptome assemblies"

    write_nonempty_gtf_list() {
      local input_dir="\$1"
      local output_list="\$2"
      local count=0

      : > "\${output_list}"
      while IFS= read -r -d '' gtf_file; do
        if [[ -s "\${gtf_file}" ]]; then
          if grep -qvE '^(#|\$)' "\${gtf_file}"; then
            printf '%s\\n' "\${gtf_file}" >> "\${output_list}"
            count=\$((count + 1))
          fi
        fi
      done < <(find "\${input_dir}" -type f -name '*.gtf' -print0 | sort -z)

      if [[ "\${count}" -eq 0 ]]; then
        return 1
      fi
      return 0
    }

    if write_nonempty_gtf_list default_gtfs default_gtfs.txt && write_nonempty_gtf_list alt_gtfs alt_gtfs.txt; then
      stringtie --merge -o merged_transcriptomes.minimap2.long_reads.default_args.gtf default_gtfs.txt
      stringtie --merge -o merged_transcriptomes.minimap2.long_reads.alt_args.gtf alt_gtfs.txt
    else
      : > merged_transcriptomes.minimap2.long_reads.default_args.gtf
      : > merged_transcriptomes.minimap2.long_reads.alt_args.gtf
    fi
    {
      printf '"%s":\n' "${task.process}"
      stringtie --version 2>&1 | awk '{ printf "  stringtie: \\"%s\\"\\n", \$0 }'
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail
    : > merged_transcriptomes.minimap2.long_reads.default_args.gtf
    : > merged_transcriptomes.minimap2.long_reads.alt_args.gtf
    printf '"%s":\n  stringtie: "stub"\n' "${task.process}" > versions.yml
    """
}
