// Quality summaries for optional ncRNA annotations produced by tRNAscan-SE
// and Infernal/Rfam. AGAT structural stats are generated per GFF3 and a
// compact custom-content TSV is passed to MultiQC.
process ncrna_annotation_qc {
  label 'process_low'

  tag "AGAT and MultiQC summaries for ncRNA annotations"
  container params.container_agat
  publishDir "${params.output_dir}/quality_report/ncrna_annotations", mode: 'copy'

  input:
    path(trna_gff3)
    path(rfam_gff3)

  output:
    path "trna_agat_stats.txt", emit: trna_agat_stats
    path "rfam_agat_stats.txt", emit: rfam_agat_stats
    path "ncrna_annotation_counts_mqc.tsv", emit: multiqc_tsv
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    run_agat_stats() {
      local input_gff="\$1"
      local output_stats="\$2"
      local label="\$3"
      if grep -qvE '^(#|\$)' "\${input_gff}"; then
        agat_sp_statistics.pl --gff "\${input_gff}" --output "\${output_stats}"
      else
        printf "No %s features available for AGAT statistics\\n" "\${label}" > "\${output_stats}"
      fi
    }

    count_features() {
      awk 'BEGIN { count = 0 } !/^#/ && NF == 9 { count += 1 } END { print count }' "\$1"
    }

    run_agat_stats "${trna_gff3}" trna_agat_stats.txt "tRNAscan-SE"
    run_agat_stats "${rfam_gff3}" rfam_agat_stats.txt "Infernal/Rfam"

    trna_count=\$(count_features "${trna_gff3}")
    rfam_count=\$(count_features "${rfam_gff3}")

    cat > ncrna_annotation_counts_mqc.tsv <<EOF
# id: titan_ncrna_annotations
# section_name: 'TITAN ncRNA annotations'
# description: 'Counts of optional tRNAscan-SE and Infernal/Rfam annotations.'
# plot_type: 'table'
Tool	GFF3 features	AGAT stats
tRNAscan-SE	\${trna_count}	trna_agat_stats.txt
Infernal/Rfam	\${rfam_count}	rfam_agat_stats.txt
EOF

    printf '"%s":\n  agat_sp_statistics: "container-pinned"\n  container: "%s"\n' \\
        "${task.process}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "Number of tRNA\\t1\\n" > trna_agat_stats.txt
    printf "Number of rRNA\\t2\\n" > rfam_agat_stats.txt
    cat > ncrna_annotation_counts_mqc.tsv <<'EOF'
# id: titan_ncrna_annotations
# section_name: 'TITAN ncRNA annotations'
# description: 'Counts of optional tRNAscan-SE and Infernal/Rfam annotations.'
# plot_type: 'table'
Tool	GFF3 features	AGAT stats
tRNAscan-SE	1	trna_agat_stats.txt
Infernal/Rfam	2	rfam_agat_stats.txt
EOF
    printf '"%s":\n  agat_sp_statistics: "stub"\n' "${task.process}" > versions.yml
    """
}
