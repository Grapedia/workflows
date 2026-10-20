// Merge split gene models with Mender (https://github.com/simrbase/mender), between the AEGIS
// annotation and the functional annotation, so every downstream step sees the merged genes.
// Only runs when params.run_mender is true (see subworkflows/aegis.nf).
process mender_split_genes {
  label 'process_aegis'
  tag "Mender: merge split gene models"

  container params.container_mender
  publishDir "${params.output_dir}/01_final_annotation/primary", mode: 'copy', saveAs: { filename ->
    if (filename in ['final_annotation.gff3', 'final_annotation_proteins_all.fasta', 'final_annotation_proteins_main.fasta']) {
      return filename
    }
    if (filename in ['mender_merges.tsv', 'mender_reverted_merges.tsv', 'mender_summary.json', 'versions.yml']) {
      return "split_gene_merging/${filename}"
    }
    if (filename.startsWith('results/titan/')) {
      return "split_gene_merging/mender/${filename.substring('results/titan/'.length())}"
    }
    if (filename == 'mender_workdir/titan/merge_candidates.txt' || filename == 'mender_workdir/titan/isoseq_validated.txt') {
      return "split_gene_merging/mender/${file(filename).name}"
    }
    return null
  }

  input:
    path(gff3, stageAs: "pre_mender_annotation.gff3")
    path(proteins_all, stageAs: "pre_mender_proteins_all.fasta")
    path(proteins_main, stageAs: "pre_mender_proteins_main.fasta")
    path(genome)
    path(reference_fastas, stageAs: "reference_proteins/*")
    path(long_read_bams, stageAs: "long_reads/*")
    val(gene_id_prefix)
    val(options)
    path(mender_script)

  output:
    path "final_annotation.gff3", emit: gff
    path "final_annotation_proteins_all.fasta", emit: proteins_all
    path "final_annotation_proteins_main.fasta", emit: proteins_main
    path "mender_merges.tsv", emit: merges_tsv
    path "mender_summary.json", emit: summary_json
    path "mender_reverted_merges.tsv", optional: true, emit: reverted_tsv
    path "results/titan/*", optional: true, emit: mender_results
    path "mender_workdir/titan/merge_candidates.txt", optional: true, emit: merge_candidates
    path "mender_workdir/titan/isoseq_validated.txt", optional: true, emit: isoseq_validated
    path "versions.yml", emit: versions

  script:
    def bams = long_read_bams instanceof List ? long_read_bams.join(' ') : long_read_bams.toString()
    def references = reference_fastas instanceof List ? reference_fastas.join(' ') : reference_fastas.toString()
    """
    set -euo pipefail
    bash ${mender_script} \\
      --gff3 ${gff3} \\
      --proteins-all ${proteins_all} \\
      --proteins-main ${proteins_main} \\
      --genome ${genome} \\
      --reference "${references}" \\
      --long-read-bams "${bams}" \\
      --threads ${task.cpus} \\
      --evalue ${options.evalue} \\
      --exclude-os "${options.exclude_os ?: ''}" \\
      --flags ${options.flags} \\
      --skip-flags ${options.skip_flags} \\
      --include-review ${options.include_review ? 'yes' : 'no'} \\
      --gene-template "${gene_id_prefix}MERGEg[GCOUNT:5]" \\
      --trans-template "${gene_id_prefix}MERGEg[GCOUNT:5]_t[TCOUNT:3]"

    printf '"%s":\\n  mender_commit: "%s"\\n  diamond: "%s"\\n  container: "%s"\\n' \\
      "${task.process}" "\$(cat /opt/mender/COMMIT 2>/dev/null || echo unknown)" \\
      "\$(diamond version | sed 's/diamond version //')" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    cp ${gff3} final_annotation.gff3
    cp ${proteins_all} final_annotation_proteins_all.fasta
    cp ${proteins_main} final_annotation_proteins_main.fasta
    printf 'merged_gene\\tchrom\\tstart\\tend\\tstrand\\tmerged_from\\n' > mender_merges.tsv
    printf '{"merged_genes": 0, "source_genes_removed": 0}\\n' > mender_summary.json
    printf '"%s":\\n  mender: "stub"\\n' "${task.process}" > versions.yml
    """
}
