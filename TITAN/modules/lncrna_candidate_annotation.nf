process lncrna_candidate_annotation {
  label 'process_transcriptome'

  tag "Preliminary lncRNA candidate annotation"
  container params.container_cpat
  publishDir "${params.output_dir}/additional_annotations/ncrna/lncrna", mode: 'copy', saveAs: { filename ->
    if (filename in [
      'lncrna_candidates.gff3',
      'lncrna_candidates.gtf',
      'lncrna_candidates.fasta',
      'lncrna_classification_summary.tsv',
      'lncrna_candidates_mqc.tsv',
      'cpat_plant.output.ORF_prob.tsv',
      'cpat_plant.output.ORF_prob.best.tsv',
      'cpat_plant.output.no_ORF.txt',
      'CPAT_run_info.log',
      'versions.yml'
    ]) {
      return filename
    }
    return null
  }

  input:
    path(genome)
    path(final_annotation_gff3)
    path(trna_gff3)
    path(rfam_gff3)
    path(star_stringtie_gtf)
    path(hisat2_stringtie_gtf)
    path(long_reads_gtf)
    path(build_lncrna_candidates)

  output:
    path "lncrna_candidates.gff3", emit: gff3
    path "lncrna_candidates.gtf", emit: gtf
    path "lncrna_candidates.fasta", emit: fasta
    path "lncrna_classification_summary.tsv", emit: summary_tsv
    path "lncrna_candidates_mqc.tsv", emit: multiqc_tsv
    path "cpat_plant.output.ORF_prob.tsv", emit: cpat_orf_prob
    path "cpat_plant.output.ORF_prob.best.tsv", emit: cpat_best
    path "cpat_plant.output.no_ORF.txt", emit: cpat_no_orf
    path "CPAT_run_info.log", emit: cpat_log
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    if [[ "${params.run_lncrna}" != "true" ]]; then
        printf "##gff-version 3\\n" > lncrna_candidates.gff3
        : > lncrna_candidates.gtf
        : > lncrna_candidates.fasta
        : > cpat_plant.output.ORF_prob.tsv
        : > cpat_plant.output.ORF_prob.best.tsv
        : > cpat_plant.output.no_ORF.txt
        : > CPAT_run_info.log
        printf "class\\tcount\\nlncRNA_candidate\\t0\\n" > lncrna_classification_summary.tsv
        cat > lncrna_candidates_mqc.tsv <<'EOF'
# id: titan_lncrna_candidates
# section_name: 'TITAN lncRNA candidates'
# description: 'Preliminary lncRNA candidates after length and coding/ncRNA overlap filters.'
# plot_type: 'table'
Metric	Count
lncRNA candidates	0
EOF
        printf '"%s":\n  lncrna_candidates: "skipped"\n  cpat_model_dir: "%s"\n  cpat_model_flavour: "%s"\n  container: "%s"\n' \\
            "${task.process}" "${params.cpat_model_dir}" "${params.cpat_model_flavour}" "${task.container}" > versions.yml
        exit 0
    fi

    if [[ "${params.lncrna_require_cpat_model}" == "true" ]]; then
        if [[ -z "${params.cpat_model_dir}" || "${params.cpat_model_dir}" == "false" ]]; then
            echo "ERROR: cpat_model_dir must be provided when run_lncrna=true and lncrna_require_cpat_model=true" >&2
            exit 1
        fi
        if [[ ! -s "${params.cpat_model_dir}/Plant_Hexamer.tsv" || ! -s "${params.cpat_model_dir}/Plant.logit.RData" ]]; then
            echo "ERROR: cpat_model_dir must contain Plant_Hexamer.tsv and Plant.logit.RData for CPAT-plant mode" >&2
            exit 1
        fi
    fi

    python3 "${build_lncrna_candidates}" \\
      --genome "${genome}" \\
      --final-annotation "${final_annotation_gff3}" \\
      --trna-gff3 "${trna_gff3}" \\
      --rfam-gff3 "${rfam_gff3}" \\
      --min-length "${params.lncrna_min_length}" \\
      --output-prefix lncrna_candidates \\
      "${star_stringtie_gtf}" "${hisat2_stringtie_gtf}" "${long_reads_gtf}"

    if [[ "${params.cpat_model_dir}" != "false" && -n "${params.cpat_model_dir}" ]] && grep -q '^>' lncrna_candidates.fasta; then
        CPAT_BIN=\$(command -v cpat.py || command -v cpat)
        "\${CPAT_BIN}" \\
          -x "${params.cpat_model_dir}/Plant_Hexamer.tsv" \\
          -d "${params.cpat_model_dir}/Plant.logit.RData" \\
          -g lncrna_candidates.fasta \\
          -o cpat_plant.output
    else
        : > cpat_plant.output.ORF_prob.tsv
        : > cpat_plant.output.ORF_prob.best.tsv
        : > cpat_plant.output.no_ORF.txt
        printf "No CPAT model or no preliminary lncRNA candidates; CPAT skipped.\\n" > CPAT_run_info.log
    fi

    python3 "${build_lncrna_candidates}" \\
      --genome "${genome}" \\
      --final-annotation "${final_annotation_gff3}" \\
      --trna-gff3 "${trna_gff3}" \\
      --rfam-gff3 "${rfam_gff3}" \\
      --min-length "${params.lncrna_min_length}" \\
      --cpat-best-tsv cpat_plant.output.ORF_prob.best.tsv \\
      --cpat-cutoff "${params.cpat_plant_cutoff}" \\
      --output-prefix lncrna_candidates \\
      "${star_stringtie_gtf}" "${hisat2_stringtie_gtf}" "${long_reads_gtf}"

    cpat_version=\$("\${CPAT_BIN:-cpat}" --version 2>&1 | head -n 1 || true)
    printf '"%s":\n  lncrna_candidates: "python-stdlib+cpat-plant"\n  cpat: "%s"\n  cpat_model_dir: "%s"\n  cpat_model_flavour: "%s"\n  cpat_cutoff: "%s"\n  output_status: "cpat_plant_filtered_candidates_not_final_vitis_annotation"\n  container: "%s"\n' \\
        "${task.process}" "\${cpat_version}" "${params.cpat_model_dir}" "${params.cpat_model_flavour}" "${params.cpat_plant_cutoff}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    cat > lncrna_candidates.gff3 <<'EOF'
##gff-version 3
chrStub	TITAN_lncRNA	gene	1	250	.	+	.	ID=lncrna_candidate_1;Name=stub_lncRNA;biotype=lncRNA_candidate
chrStub	TITAN_lncRNA	lnc_RNA	1	250	.	+	.	ID=lncrna_candidate_1.t1;Parent=lncrna_candidate_1;Name=stub_lncRNA;biotype=lncRNA_candidate
chrStub	TITAN_lncRNA	exon	1	250	.	+	.	ID=lncrna_candidate_1.t1.exon1;Parent=lncrna_candidate_1.t1
EOF
    printf 'chrStub\\tTITAN_lncRNA\\texon\\t1\\t250\\t.\\t+\\t.\\tgene_id "lncrna_candidate_1"; transcript_id "lncrna_candidate_1.t1";\\n' > lncrna_candidates.gtf
    printf ">lncrna_candidate_1.t1 length=250\\n%s\\n" "ACGTACGTACGT" > lncrna_candidates.fasta
    printf "seq_ID\\tcoding_prob\\nlncrna_candidate_1.t1\\t0.12\\n" > cpat_plant.output.ORF_prob.best.tsv
    cp cpat_plant.output.ORF_prob.best.tsv cpat_plant.output.ORF_prob.tsv
    : > cpat_plant.output.no_ORF.txt
    printf "CPAT stub\\n" > CPAT_run_info.log
    printf "class\\tcount\\nlncRNA_candidate\\t1\\nexcluded_cpat_coding\\t0\\n" > lncrna_classification_summary.tsv
    cat > lncrna_candidates_mqc.tsv <<'EOF'
# id: titan_lncrna_candidates
# section_name: 'TITAN lncRNA candidates'
# description: 'Preliminary lncRNA candidates after length and coding/ncRNA overlap filters.'
# plot_type: 'table'
Metric	Count
lncRNA candidates	1
EOF
    printf '"%s":\n  lncrna_candidates: "stub"\n' "${task.process}" > versions.yml
    """
}
