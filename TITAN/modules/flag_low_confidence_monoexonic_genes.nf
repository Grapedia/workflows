process flag_low_confidence_monoexonic_genes {
  label 'process_low'
  tag "Classify single-exon genes by supporting evidence"

  container params.container_python
  publishDir "${params.output_dir}/quality_report/monoexonic_gene_confidence", mode: 'copy'

  input:
    path(annotation)
    path(expression_json)
    path(te_gff3)
    path(interproscan_tsv)
    path(eggnog_tsv)
    path(diamond2go_tsv)
    path(liftoff_correspondence_tsv)
    path(busco_full_table)
    path(proteins_main_fasta)
    path(script)

  output:
    path "monoexonic_gene_confidence.tsv", emit: tsv_report
    path "monoexonic_gene_confidence_summary.json", emit: json_summary
    path "final_annotation.high_confidence.gff3", emit: filtered_gff3
    path "final_annotation_proteins_main.high_confidence.fasta", emit: filtered_proteins_main
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    python3 ${script} \\
      --gff3 ${annotation} \\
      --expression-json ${expression_json} \\
      --te-gff3 ${te_gff3} \\
      --interproscan-tsv ${interproscan_tsv} \\
      --eggnog-tsv ${eggnog_tsv} \\
      --diamond2go-tsv ${diamond2go_tsv} \\
      --liftoff-correspondence-tsv ${liftoff_correspondence_tsv} \\
      --busco-full-table ${busco_full_table} \\
      --proteins-main-fasta ${proteins_main_fasta} \\
      --tsv-report monoexonic_gene_confidence.tsv \\
      --json-summary monoexonic_gene_confidence_summary.json \\
      --filtered-gff3 final_annotation.high_confidence.gff3 \\
      --filtered-proteins-main final_annotation_proteins_main.high_confidence.fasta
    script_sha256=\$(sha256sum ${script} | awk '{print \$1}')
    printf '"%s":\n  container: "%s"\n  script: "%s"\n  script_sha256: "%s"\n' \\
      "${task.process}" "${task.container}" "${script}" "\${script_sha256}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    python3 ${script} \\
      --gff3 ${annotation} \\
      --expression-json ${expression_json} \\
      --te-gff3 ${te_gff3} \\
      --interproscan-tsv ${interproscan_tsv} \\
      --eggnog-tsv ${eggnog_tsv} \\
      --diamond2go-tsv ${diamond2go_tsv} \\
      --liftoff-correspondence-tsv ${liftoff_correspondence_tsv} \\
      --busco-full-table ${busco_full_table} \\
      --proteins-main-fasta ${proteins_main_fasta} \\
      --tsv-report monoexonic_gene_confidence.tsv \\
      --json-summary monoexonic_gene_confidence_summary.json \\
      --filtered-gff3 final_annotation.high_confidence.gff3 \\
      --filtered-proteins-main final_annotation_proteins_main.high_confidence.fasta
    script_sha256=\$(sha256sum ${script} | awk '{print \$1}')
    printf '"%s":\n  container: "%s"\n  script: "%s"\n  script_sha256: "%s"\n' \\
      "${task.process}" "${task.container}" "${script}" "\${script_sha256}" > versions.yml
    """
}
