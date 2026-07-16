process aegis_short_reads {
  label 'process_aegis'
    tag "AEGIS merge of short-read evidence"

    container "${params.aegis_container}"
    publishDir "${params.output_dir}/aegis_outputs", mode: 'copy'
    input:
        path(edta_masked_genome)
        path(augustus_gff)
        path(genemark_gtf)
        path(liftoff_annotations)
        path(egapx_gff3)
        path(stranded_default_args)
        path(stranded_alt_args)
        path(gffcompare_stranded)
        path(gffcompare_unstranded)
        path(unstranded_default_args)
        path(unstranded_alt_args)
        path(aegis_merge_script)

    output:
        path "final_annotation.gff3", emit: aegis_gff
        path "final_annotation_proteins_all.fasta", emit: aegis_proteins_all
        path "final_annotation_proteins_main.fasta", emit: aegis_proteins_main
    path "versions.yml", emit: versions

    script:
        """
        set -euo pipefail
        export NXF_TASK_PROCESS="${task.process}"
        bash ${aegis_merge_script} \\
          short_reads \\
          ${edta_masked_genome} \\
          ${params.aegis_version} \\
          ${params.aegis_container} \\
          final_annotation \\
          ${liftoff_annotations} \\
          ${augustus_gff} \\
          ${genemark_gtf} \\
          ${egapx_gff3} \\
          ${stranded_default_args} \\
          ${stranded_alt_args} \\
          ${gffcompare_stranded} \\
          ${gffcompare_unstranded} \\
          ${unstranded_default_args} \\
          ${unstranded_alt_args}
        """

    stub:
        """
        seqid=\$(awk '/^>/ {print substr(\$1, 2); exit}' ${edta_masked_genome})
        length=\$(awk 'BEGIN{n=0} !/^>/ {n += length(\$0)} END{print n}' ${edta_masked_genome})
        end=\$(( length < 9 ? length : 9 ))
        printf "##gff-version 3\\n%s\\tAegis\\tgene\\t1\\t%s\\t.\\t+\\t.\\tID=aegis_stub_gene\\n%s\\tAegis\\tmRNA\\t1\\t%s\\t.\\t+\\t.\\tID=aegis_stub_tx;Parent=aegis_stub_gene\\n%s\\tAegis\\tCDS\\t1\\t%s\\t.\\t+\\t0\\tID=aegis_stub_cds;Parent=aegis_stub_tx\\n" "\$seqid" "\$end" "\$seqid" "\$end" "\$seqid" "\$end" > final_annotation.gff3
        printf ">aegis_stub_protein\\nM\\n" > final_annotation_proteins_all.fasta
        printf ">aegis_stub_protein\\nM\\n" > final_annotation_proteins_main.fasta
        printf '"%s":\\n  aegis: "%s"\\n  aegis_container: "%s"\\n' \\
          "${task.process}" "${params.aegis_version}" "${params.aegis_container}" > versions.yml
        """
}
