// Aggregates run-wide QC (fastp trimming reports across every RNA-seq
// sample) with final-annotation quality signals (BUSCO completeness, AGAT
// structural stats, optional ncRNA/SQANTI3 annotation QC and the structural
// validation report) into one HTML report.
process multiqc_report {
  label 'process_low'

  tag "MultiQC full-run and final annotation quality report"
  container params.container_multiqc
  publishDir "${params.output_dir}/quality_report", mode: 'copy'

  input:
    path(fastp_json_reports, stageAs: "fastp_reports/*")
    path(busco_short_summary)
    path(omark_mqc_tsv)
    path(agat_stats_txt)
    path(ncrna_qc_reports, stageAs: "ncrna_qc/*")
    path(lncrna_qc_tsv)
    path(sqanti3_qc_tsv)
    path(expression_support_mqc_tsv)
    path(final_annotation_sources_qc_tsv)
    path(validation_json)

  output:
    path "titan_multiqc_report.html", emit: html_report
    path "titan_multiqc_report_data", emit: report_data
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    mkdir -p mqc_input
    cp fastp_reports/*.fastp.json mqc_input/ 2>/dev/null || true
    cp ${busco_short_summary} mqc_input/ 2>/dev/null || true
    cp ${omark_mqc_tsv} mqc_input/ 2>/dev/null || true
    cp ${agat_stats_txt} mqc_input/ 2>/dev/null || true
    cp ncrna_qc/* mqc_input/ 2>/dev/null || true
    cp ${lncrna_qc_tsv} mqc_input/ 2>/dev/null || true
    cp ${sqanti3_qc_tsv} mqc_input/ 2>/dev/null || true
    cp ${expression_support_mqc_tsv} mqc_input/ 2>/dev/null || true
    cp ${final_annotation_sources_qc_tsv} mqc_input/ 2>/dev/null || true
    cp ${validation_json} mqc_input/ 2>/dev/null || true

    multiqc mqc_input --filename titan_multiqc_report.html --force

    multiqc --version | sed 's/^/  multiqc: "/; s/\$/"/' | {
      printf '"%s":\\n' "${task.process}"
      cat
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "<html>stub</html>\\n" > titan_multiqc_report.html
    mkdir -p titan_multiqc_report_data
    printf '"%s":\n  multiqc: "stub"\n' "${task.process}" > versions.yml
    """
}
