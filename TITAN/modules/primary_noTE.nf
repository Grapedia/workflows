// Protein-coding, TE-cleaned variant of the final annotation (01_final_annotation/primary_noTE).
// Three steps: select the TE-derived candidate genes, test them for a non-Vitis homolog with
// DIAMOND (candidates only, a few thousand proteins), then write the filtered GFF3 / proteins /
// functional-annotation tables. The last step also removes the genes carrying a TE functional
// annotation listed in params.noTE_te_function_terms (Liftoff genes kept).
// See scripts/filter_te_ncrna_final_annotation.py.

process select_te_candidates {
  label 'process_low'
  tag "Select TE-derived gene candidates (CDS >= ${params.noTE_te_frac} in EDTA TE)"

  container params.container_python

  input:
    path(annotation)
    path(proteins_main_fasta)
    path(te_gff3)
    path(reference_gff3)
    path(liftoff_correspondence_tsv)
    path(busco_full_table)
    path(diamond2go_main_tsv)
    path(eggnog_main_annotations)
    path(interproscan_main_tsv)
    path(script)

  output:
    path "te_candidates.tsv", emit: candidates_tsv
    path "te_candidates.fasta", emit: candidates_fasta
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    python3 ${script} select \\
      --gff3 ${annotation} \\
      --proteins-main ${proteins_main_fasta} \\
      --te-gff3 ${te_gff3} \\
      --reference-gff3 ${reference_gff3} \\
      --liftoff-correspondence-tsv ${liftoff_correspondence_tsv} \\
      --busco-full-table ${busco_full_table} \\
      --diamond2go-main ${diamond2go_main_tsv} \\
      --eggnog-main ${eggnog_main_annotations} \\
      --interproscan-main ${interproscan_main_tsv} \\
      --te-frac ${params.noTE_te_frac} \\
      --out-dir .
    script_sha256=\$(sha256sum ${script} | awk '{print \$1}')
    printf '"%s":\\n  container: "%s"\\n  script: "%s"\\n  script_sha256: "%s"\\n' \\
      "${task.process}" "${task.container}" "${script}" "\${script_sha256}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    python3 ${script} select \\
      --gff3 ${annotation} \\
      --proteins-main ${proteins_main_fasta} \\
      --te-gff3 ${te_gff3} \\
      --reference-gff3 ${reference_gff3} \\
      --liftoff-correspondence-tsv ${liftoff_correspondence_tsv} \\
      --busco-full-table ${busco_full_table} \\
      --diamond2go-main ${diamond2go_main_tsv} \\
      --eggnog-main ${eggnog_main_annotations} \\
      --interproscan-main ${interproscan_main_tsv} \\
      --te-frac ${params.noTE_te_frac} \\
      --out-dir .
    printf '"%s":\\n  container: "not_recorded"\\n' "${task.process}" > versions.yml
    """
}

process diamond_te_candidates {
  label 'process_aegis'
  tag "DIAMOND of the TE-derived candidates against non-Vitis protein evidence"

  container params.container_diamond2go

  input:
    path(candidates_fasta)
    path(protein_evidence_fastas)

  output:
    path "te_candidates_vs_nonvitis.tsv", emit: hits_tsv
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    # Non-Vitis reference: UniProt-style records (header carries 'OS=') of the protein samplesheet
    # FASTAs, minus 'OS=Vitis'. A hit in another genus is evidence independent of any Vitis annotation.
    cat ${protein_evidence_fastas} \\
      | awk '/^>/{keep = (\$0 ~ /OS=/ && \$0 !~ /OS=Vitis /)} keep' \\
      | sed '/^>/! s/[*-]//g' > nonvitis_reference.fasta
    if grep -q '>' nonvitis_reference.fasta && grep -q '>' ${candidates_fasta}; then
      diamond makedb --in nonvitis_reference.fasta -d nonvitis_reference --threads ${task.cpus}
      diamond blastp -q ${candidates_fasta} -d nonvitis_reference -o te_candidates_vs_nonvitis.tsv \\
        --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore \\
        --threads ${task.cpus} --sensitive -e 1e-5 -k 3
    else
      echo "WARNING: no UniProt-style (OS=) non-Vitis protein records or no candidates; no homology rescue" >&2
      : > te_candidates_vs_nonvitis.tsv
    fi
    printf '"%s":\\n  container: "%s"\\n  diamond: "%s"\\n' \\
      "${task.process}" "${task.container}" "\$(diamond version | sed 's/diamond version //')" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    : > te_candidates_vs_nonvitis.tsv
    printf '"%s":\\n  diamond: "stub"\\n' "${task.process}" > versions.yml
    """
}

process finalize_primary_noTE {
  label 'process_low'
  tag "Write primary_noTE (TE-derived genes, TE-function genes and ncRNA removed)"

  container params.container_python
  publishDir "${params.output_dir}/01_final_annotation/primary_noTE", mode: 'copy'

  input:
    path(annotation)
    path(proteins_all_fasta)
    path(proteins_main_fasta)
    path(candidates_tsv)
    path(diamond_hits_tsv)
    path(diamond2go_all_tsv)
    path(diamond2go_main_tsv)
    path(eggnog_all_annotations)
    path(eggnog_main_annotations)
    path(interproscan_all_tsv)
    path(interproscan_main_tsv)
    path(te_function_terms)
    path(liftoff_correspondence_tsv)
    path(script)

  output:
    path "final_annotation.noTE.gff3", emit: gff3
    path "final_annotation_proteins_all.noTE.fasta", emit: proteins_all
    path "final_annotation_proteins_main.noTE.fasta", emit: proteins_main
    path "functional_annotation", emit: functional_annotation
    path "genes_removed.tsv", emit: removed_tsv
    path "primary_noTE_summary.json", emit: summary_json
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    python3 ${script} finalize \\
      --gff3 ${annotation} \\
      --proteins-all ${proteins_all_fasta} \\
      --proteins-main ${proteins_main_fasta} \\
      --candidates-tsv ${candidates_tsv} \\
      --diamond-hits ${diamond_hits_tsv} \\
      --diamond2go-all ${diamond2go_all_tsv} \\
      --diamond2go-main ${diamond2go_main_tsv} \\
      --eggnog-all ${eggnog_all_annotations} \\
      --eggnog-main ${eggnog_main_annotations} \\
      --interproscan-all ${interproscan_all_tsv} \\
      --interproscan-main ${interproscan_main_tsv} \\
      --te-function-terms ${te_function_terms} \\
      --liftoff-correspondence-tsv ${liftoff_correspondence_tsv} \\
      --te-frac ${params.noTE_te_frac} \\
      --out-dir .
    script_sha256=\$(sha256sum ${script} | awk '{print \$1}')
    printf '"%s":\\n  container: "%s"\\n  script: "%s"\\n  script_sha256: "%s"\\n' \\
      "${task.process}" "${task.container}" "${script}" "\${script_sha256}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    python3 ${script} finalize \\
      --gff3 ${annotation} \\
      --proteins-all ${proteins_all_fasta} \\
      --proteins-main ${proteins_main_fasta} \\
      --candidates-tsv ${candidates_tsv} \\
      --diamond-hits ${diamond_hits_tsv} \\
      --diamond2go-all ${diamond2go_all_tsv} \\
      --diamond2go-main ${diamond2go_main_tsv} \\
      --eggnog-all ${eggnog_all_annotations} \\
      --eggnog-main ${eggnog_main_annotations} \\
      --interproscan-all ${interproscan_all_tsv} \\
      --interproscan-main ${interproscan_main_tsv} \\
      --te-function-terms ${te_function_terms} \\
      --liftoff-correspondence-tsv ${liftoff_correspondence_tsv} \\
      --te-frac ${params.noTE_te_frac} \\
      --out-dir .
    printf '"%s":\\n  container: "not_recorded"\\n' "${task.process}" > versions.yml
    """
}

// Side-by-side quality table of the final-annotation candidates (primary, primary_noTE,
// high_confidence_monoexonic): gene counts, TE content, function, expression, BUSCO, OMArk.
process compare_candidates {
  label 'process_low'
  tag "Compare the final annotation candidates (primary / primary_noTE / high_confidence_monoexonic)"

  container params.container_python
  publishDir "${params.output_dir}/01_final_annotation/quality_report", mode: 'copy'

  input:
    path(primary_gff3, stageAs: "primary/*")
    path(primary_busco, stageAs: "primary/*")
    path(primary_omark, stageAs: "primary/*")
    path(noTE_gff3, stageAs: "primary_noTE/*")
    path(noTE_busco, stageAs: "primary_noTE/*")
    path(noTE_omark, stageAs: "primary_noTE/*")
    path(monoexonic_gff3, stageAs: "high_confidence_monoexonic/*")
    path(monoexonic_busco, stageAs: "high_confidence_monoexonic/*")
    path(te_gff3)
    path(expression_json)
    path(diamond2go_main_tsv)
    path(eggnog_main_annotations)
    path(interproscan_main_tsv)
    path(script)

  output:
    path "candidates_comparison.tsv", emit: comparison_tsv
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    python3 ${script} compare \\
      --candidate primary ${primary_gff3} ${primary_busco} ${primary_omark} \\
      --candidate primary_noTE ${noTE_gff3} ${noTE_busco} ${noTE_omark} \\
      --candidate high_confidence_monoexonic ${monoexonic_gff3} ${monoexonic_busco} NA \\
      --te-gff3 ${te_gff3} \\
      --expression-json ${expression_json} \\
      --diamond2go-main ${diamond2go_main_tsv} \\
      --eggnog-main ${eggnog_main_annotations} \\
      --interproscan-main ${interproscan_main_tsv} \\
      --out candidates_comparison.tsv
    script_sha256=\$(sha256sum ${script} | awk '{print \$1}')
    printf '"%s":\\n  container: "%s"\\n  script: "%s"\\n  script_sha256: "%s"\\n' \\
      "${task.process}" "${task.container}" "${script}" "\${script_sha256}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    python3 ${script} compare \\
      --candidate primary ${primary_gff3} ${primary_busco} ${primary_omark} \\
      --candidate primary_noTE ${noTE_gff3} ${noTE_busco} ${noTE_omark} \\
      --candidate high_confidence_monoexonic ${monoexonic_gff3} ${monoexonic_busco} NA \\
      --te-gff3 ${te_gff3} \\
      --expression-json ${expression_json} \\
      --diamond2go-main ${diamond2go_main_tsv} \\
      --eggnog-main ${eggnog_main_annotations} \\
      --interproscan-main ${interproscan_main_tsv} \\
      --out candidates_comparison.tsv
    printf '"%s":\\n  container: "not_recorded"\\n' "${task.process}" > versions.yml
    """
}
