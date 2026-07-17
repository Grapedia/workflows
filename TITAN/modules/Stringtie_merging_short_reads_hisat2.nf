// Transcriptomes merging with StringTie
process Stringtie_merging_short_reads_hisat2 {
  label 'process_merge'

  tag "HISAT2/StringTie merge: stranded and optional unstranded short-read GTFs"
  container params.container_stringtie
  publishDir "${params.output_dir}/tmp", mode: 'copy', enabled: params.publish_intermediates
  publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename ->
    if (filename == 'merged_transcriptomes.hisat2.short_reads.default_args.stranded.gtf') {
      return 'merged_hisat2_stringtie_stranded_default.gtf'
    }
    if (filename == 'merged_transcriptomes.hisat2.short_reads.alt_args.stranded.gtf') {
      return 'merged_hisat2_stringtie_stranded_alt.gtf'
    }
    if (filename == 'merged_transcriptomes.hisat2.short_reads.default_args.unstranded.gtf') {
      return 'merged_hisat2_stringtie_unstranded_default.gtf'
    }
    if (filename == 'merged_transcriptomes.hisat2.short_reads.alt_args.unstranded.gtf') {
      return 'merged_hisat2_stringtie_unstranded_alt.gtf'
    }
    return null
  }
  input:
    path(stranded_default_gtfs, stageAs: "stranded_default_gtfs/*")
    path(stranded_alt_gtfs, stageAs: "stranded_alt_gtfs/*")
    path(unstranded_default_gtfs, stageAs: "unstranded_default_gtfs/*")
    path(unstranded_alt_gtfs, stageAs: "unstranded_alt_gtfs/*")

  output:
    path "merged_transcriptomes.hisat2.short_reads.default_args.stranded.gtf", emit: default_args_stranded
    path "merged_transcriptomes.hisat2.short_reads.alt_args.stranded.gtf", emit: alt_args_stranded
    path "merged_transcriptomes.hisat2.short_reads.default_args.unstranded.gtf", emit: default_args_unstranded
    path "merged_transcriptomes.hisat2.short_reads.alt_args.unstranded.gtf", emit: alt_args_unstranded
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running StringTie merging on HISAT2/StringTie transcriptomes - separating stranded and unstranded samples."

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

    write_nonempty_gtf_list stranded_default_gtfs stranded_default_gtfs.txt true
    write_nonempty_gtf_list stranded_alt_gtfs stranded_alt_gtfs.txt true
    stringtie --merge -o merged_transcriptomes.hisat2.short_reads.default_args.stranded.gtf stranded_default_gtfs.txt
    stringtie --merge -o merged_transcriptomes.hisat2.short_reads.alt_args.stranded.gtf stranded_alt_gtfs.txt

    if write_nonempty_gtf_list unstranded_default_gtfs unstranded_default_gtfs.txt false && write_nonempty_gtf_list unstranded_alt_gtfs unstranded_alt_gtfs.txt false; then
        echo "[\$DATE] Running StringTie merging on HISAT2/StringTie transcriptomes - unstranded samples detected."
        stringtie --merge -o merged_transcriptomes.hisat2.short_reads.default_args.unstranded.gtf unstranded_default_gtfs.txt
        stringtie --merge -o merged_transcriptomes.hisat2.short_reads.alt_args.unstranded.gtf unstranded_alt_gtfs.txt
    else
        : > merged_transcriptomes.hisat2.short_reads.default_args.unstranded.gtf
        : > merged_transcriptomes.hisat2.short_reads.alt_args.unstranded.gtf
    fi
    {
      printf '"%s":\n' "${task.process}"
      stringtie --version 2>&1 | awk '{ printf "  stringtie: \\"%s\\"\\n", \$0 }'
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail
    : > merged_transcriptomes.hisat2.short_reads.default_args.stranded.gtf
    : > merged_transcriptomes.hisat2.short_reads.alt_args.stranded.gtf
    : > merged_transcriptomes.hisat2.short_reads.default_args.unstranded.gtf
    : > merged_transcriptomes.hisat2.short_reads.alt_args.unstranded.gtf
    printf '"%s":\n  stringtie: "stub"\n' "${task.process}" > versions.yml
    """
}
