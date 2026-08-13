process aegis_liftoff_gene_ids {
  label 'process_aegis'
  tag "Carry over liftoff gene IDs onto the final annotation"

  container "${params.container_aegis}"
  publishDir "${params.output_dir}/aegis_outputs", mode: 'copy', saveAs: { filename ->
    if (filename in ['final_annotation.gff3', 'final_annotation_proteins_all.fasta', 'final_annotation_proteins_main.fasta', 'liftoff_gene_id_correspondence.tsv', 'versions.yml']) {
      return filename
    }
    return null
  }
  publishDir "${params.output_dir}/intermediate_files/aegis", mode: 'copy', enabled: params.publish_intermediates, saveAs: { filename ->
    if (filename == 'overlap_out') {
      return filename
    }
    return null
  }

  input:
    path(vitvi_gff3, stageAs: "vitvi_annotation.gff3")
    path(vitvi_proteins_all, stageAs: "vitvi_proteins_all.fasta")
    path(vitvi_proteins_main, stageAs: "vitvi_proteins_main.fasta")
    path(liftoff_gff3)
    path(previous_annotations)
    path(apply_liftoff_ids_script)

  output:
    path "final_annotation.gff3", emit: aegis_gff
    path "final_annotation_proteins_all.fasta", emit: aegis_proteins_all
    path "final_annotation_proteins_main.fasta", emit: aegis_proteins_main
    path "liftoff_gene_id_correspondence.tsv", emit: correspondence_tsv
    path "overlap_out", optional: true, emit: overlap_dir
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    mkdir -p overlap_out
    /opt/conda/envs/bio_env/bin/python -m aegis overlap \\
      -a final,liftoff \\
      --original-annotation-files "NA,${previous_annotations}" \\
      -r final \\
      -d overlap_out \\
      ${vitvi_gff3} ${liftoff_gff3}

    overlap_csv=\$(find overlap_out -name '*_overlaps_t*.csv' -print -quit)
    if [[ -z "\${overlap_csv:-}" ]]; then
      echo "No aegis overlap correspondence produced; keeping Vitvi IDs unchanged" >&2
      cp ${vitvi_gff3} final_annotation.gff3
      cp ${vitvi_proteins_all} final_annotation_proteins_all.fasta
      cp ${vitvi_proteins_main} final_annotation_proteins_main.fasta
      printf "original_gene_id\\tfinal_gene_id\\tdecision\\toverlap_score\\tmin_gene_percent\\tmin_cds_percent\\n" > liftoff_gene_id_correspondence.tsv
    else
      python3 ${apply_liftoff_ids_script} \\
        --overlap-csv "\$overlap_csv" \\
        --gff3 ${vitvi_gff3} \\
        --proteins-all ${vitvi_proteins_all} \\
        --proteins-main ${vitvi_proteins_main} \\
        --new-origin final \\
        --old-origin liftoff \\
        -o final_annotation.gff3 \\
        --output-proteins-all final_annotation_proteins_all.fasta \\
        --output-proteins-main final_annotation_proteins_main.fasta \\
        --audit-tsv liftoff_gene_id_correspondence.tsv
    fi

    test -s final_annotation.gff3
    test -s final_annotation_proteins_all.fasta
    test -s final_annotation_proteins_main.fasta

    printf '"%s":\\n  aegis_overlap: "used"\\n  container: "%s"\\n' "${task.process}" "${params.container_aegis}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    cp ${vitvi_gff3} final_annotation.gff3
    cp ${vitvi_proteins_all} final_annotation_proteins_all.fasta
    cp ${vitvi_proteins_main} final_annotation_proteins_main.fasta
    printf "original_gene_id\\tfinal_gene_id\\tdecision\\toverlap_score\\tmin_gene_percent\\tmin_cds_percent\\n" > liftoff_gene_id_correspondence.tsv
    printf '"%s":\\n  aegis_overlap: "stub"\\n' "${task.process}" > versions.yml
    """
}
