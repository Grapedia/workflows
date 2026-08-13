process aegis_merge {
  label 'process_aegis'
  tag "AEGIS rename/tidy of mikado_pick annotation"

  // Output here still carries systematic Vitvi IDs for every gene; it is an
  // intermediate product until aegis_liftoff_gene_ids carries over matching
  // old liftoff IDs, which is what actually gets published as the final
  // annotation under aegis_outputs.
  container "${params.container_aegis}"
  publishDir "${params.output_dir}/intermediate_files/aegis", mode: 'copy', enabled: params.publish_intermediates, saveAs: { filename ->
    if (filename in ['final_annotation.gff3', 'final_annotation_proteins_all.fasta', 'final_annotation_proteins_main.fasta', 'aegis_rename', 'aegis_tidy', 'aegis_proteins_all', 'aegis_proteins_main', 'aegis_inputs.tsv', 'aegis_merge.log', 'versions.yml']) {
      return "vitvi_only_${filename}"
    }
    return null
  }

  input:
    path(edta_masked_genome)
    path(mikado_gff3)
    path(aegis_merge_script)

  output:
    path "final_annotation.gff3", emit: aegis_gff
    path "final_annotation_proteins_all.fasta", emit: aegis_proteins_all
    path "final_annotation_proteins_main.fasta", emit: aegis_proteins_main
    path "aegis_inputs.tsv", emit: input_manifest
    path "aegis_merge.log", emit: debug_log
    path "aegis_rename", optional: true, emit: aegis_rename_dir
    path "aegis_tidy", optional: true, emit: aegis_tidy_dir
    path "aegis_proteins_all", optional: true, emit: aegis_proteins_all_dir
    path "aegis_proteins_main", optional: true, emit: aegis_proteins_main_dir
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    export NXF_TASK_PROCESS="${task.process}"
    bash ${aegis_merge_script} \\
      ${edta_masked_genome} \\
      ${params.aegis_version} \\
      ${params.container_aegis} \\
      final_annotation \\
      ${mikado_gff3} \\
      ${params.aegis_gene_id_prefix} 2>&1 | tee aegis_merge.log
    """

  stub:
    """
    set -euo pipefail
    seqid=\$(awk '/^>/ {print substr(\$1, 2); exit}' ${edta_masked_genome})
    length=\$(awk 'BEGIN{n=0} !/^>/ {n += length(\$0)} END{print n}' ${edta_masked_genome})
    end=\$(( length < 9 ? length : 9 ))
    printf "##gff-version 3\\n%s\\tAegis\\tgene\\t1\\t%s\\t.\\t+\\t.\\tID=aegis_stub_gene\\n%s\\tAegis\\tmRNA\\t1\\t%s\\t.\\t+\\t.\\tID=aegis_stub_tx;Parent=aegis_stub_gene\\n%s\\tAegis\\tCDS\\t1\\t%s\\t.\\t+\\t0\\tID=aegis_stub_cds;Parent=aegis_stub_tx\\n" "\$seqid" "\$end" "\$seqid" "\$end" "\$seqid" "\$end" > final_annotation.gff3
    printf ">aegis_stub_protein\\nM\\n" > final_annotation_proteins_all.fasta
    printf ">aegis_stub_protein\\nM\\n" > final_annotation_proteins_main.fasta
    printf "name\\tpath\\trequired\\tincluded\\tsize_bytes\\n" > aegis_inputs.tsv
    printf "Stub AEGIS finalize\\n" > aegis_merge.log
    printf '"%s":\\n  aegis: "%s"\\n  aegis_container: "%s"\\n' \\
      "${task.process}" "${params.aegis_version}" "${params.container_aegis}" > versions.yml
    """
}
