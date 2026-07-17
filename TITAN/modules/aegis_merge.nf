process aegis_merge {
  label 'process_aegis'
  tag "AEGIS merge of ${aegis_mode} evidence"

  container "${params.container_aegis}"
  publishDir "${params.output_dir}/aegis_outputs", mode: 'copy', saveAs: { filename ->
    if (filename in ['final_annotation.gff3', 'final_annotation_proteins_all.fasta', 'final_annotation_proteins_main.fasta', 'versions.yml']) {
      return filename
    }
    return null
  }
  publishDir "${params.output_dir}/intermediate_files/aegis", mode: 'copy', enabled: params.publish_intermediates, saveAs: { filename ->
    if (filename in ['aegis_merge', 'aegis_rename', 'aegis_tidy', 'aegis_proteins_all', 'aegis_proteins_main', 'aegis_inputs.tsv', 'aegis_merge.log']) {
      return filename
    }
    return null
  }

  input:
    val(aegis_mode)
    path(edta_masked_genome)
    path(augustus_gff)
    path(genemark_gtf)
    path(liftoff_annotations)
    path(egapx_gff3)
    path(long_reads_default_args)
    path(long_reads_alt_args)
    path(stranded_default_args)
    path(stranded_alt_args)
    path(gffcompare_stranded)
    path(gffcompare_unstranded)
    path(unstranded_default_args)
    path(unstranded_alt_args)
    path(helixer_gff3)
    path(aegis_merge_script)

  output:
    path "final_annotation.gff3", emit: aegis_gff
    path "final_annotation_proteins_all.fasta", emit: aegis_proteins_all
    path "final_annotation_proteins_main.fasta", emit: aegis_proteins_main
    path "aegis_inputs.tsv", emit: input_manifest
    path "aegis_merge.log", emit: debug_log
    path "aegis_merge", optional: true, emit: aegis_merge_dir
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
      ${aegis_mode} \\
      ${edta_masked_genome} \\
      ${params.aegis_version} \\
      ${params.container_aegis} \\
      final_annotation \\
      ${liftoff_annotations} \\
      ${augustus_gff} \\
      ${genemark_gtf} \\
      ${egapx_gff3} \\
      ${stranded_default_args} \\
      ${stranded_alt_args} \\
      ${gffcompare_stranded} \\
      ${long_reads_default_args} \\
      ${long_reads_alt_args} \\
      ${gffcompare_unstranded} \\
      ${unstranded_default_args} \\
      ${unstranded_alt_args} \\
      ${helixer_gff3} \\
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
    printf "Stub AEGIS merge\\n" > aegis_merge.log
    printf '"%s":\\n  aegis: "%s"\\n  aegis_container: "%s"\\n' \\
      "${task.process}" "${params.aegis_version}" "${params.container_aegis}" > versions.yml
    """
}
