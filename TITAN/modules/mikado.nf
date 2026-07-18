process mikado_prepare {
  label 'process_transcriptome'
  tag "Mikado configure and prepare"
  container params.container_mikado
  publishDir "${params.output_dir}/final_annotations/mikado/intermediate", mode: 'copy', enabled: params.publish_intermediates, saveAs: { filename ->
    if (filename in ['mikado_configuration.yaml', 'mikado_prepared.fasta', 'mikado_prepared.gtf', 'transcript_inputs.tsv', 'versions.yml']) {
      return filename
    }
    return null
  }

  input:
    path(genome)
    path(braker_augustus_gff3)
    path(braker_genemark_gtf)
    path(liftoff_gff3)
    path(egapx_gff3)
    path(star_stringtie_default_stranded)
    path(star_stringtie_alt_stranded)
    path(star_psiclass_stranded)
    path(star_psiclass_unstranded)
    path(star_stringtie_default_unstranded)
    path(star_stringtie_alt_unstranded)
    path(hisat2_stringtie_default_stranded)
    path(hisat2_stringtie_alt_stranded)
    path(hisat2_stringtie_default_unstranded)
    path(hisat2_stringtie_alt_unstranded)
    path(long_reads_default)
    path(long_reads_alt)
    path(flair_isoforms_gtf)
    path(helixer_gff3)
    path(make_mikado_list)

  output:
    path "mikado_configuration.yaml", emit: config
    path "mikado_prepared.fasta", emit: fasta
    path "mikado_prepared.gtf", emit: gtf
    path "transcript_inputs.tsv", emit: input_list
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    if [[ "${params.run_mikado}" != "true" ]]; then
      printf "# Mikado skipped; set --run_mikado true to enable.\\n" > mikado_configuration.yaml
      : > mikado_prepared.fasta
      : > mikado_prepared.gtf
      printf "path\\tlabel\\tstranded\\tscore\\tis_reference\\n" > transcript_inputs.tsv
      printf '"%s":\n  mikado_prepare: "skipped"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
      exit 0
    fi

    # Keep only valid GFF3/comment lines before Mikado's strict parser.
    awk -F'\\t' 'NF == 9 || /^#/' "${liftoff_gff3}" > liftoff_sanitized.gff3

    python3 "${make_mikado_list}" \\
      --source "liftoff_sanitized.gff3:liftoff:False:10:True" \\
      --source "${egapx_gff3}:egapx:False:9:True" \\
      --source "${braker_augustus_gff3}:braker_augustus:False:8:False" \\
      --source "${braker_genemark_gtf}:braker_genemark:False:7:False" \\
      --source "${star_stringtie_default_stranded}:star_stringtie_default_stranded:True:5:False" \\
      --source "${star_stringtie_alt_stranded}:star_stringtie_alt_stranded:True:4:False" \\
      --source "${star_psiclass_stranded}:star_psiclass_stranded:True:5:False" \\
      --source "${star_psiclass_unstranded}:star_psiclass_unstranded:False:3:False" \\
      --source "${star_stringtie_default_unstranded}:star_stringtie_default_unstranded:False:3:False" \\
      --source "${star_stringtie_alt_unstranded}:star_stringtie_alt_unstranded:False:2:False" \\
      --source "${hisat2_stringtie_default_stranded}:hisat2_stringtie_default_stranded:True:5:False" \\
      --source "${hisat2_stringtie_alt_stranded}:hisat2_stringtie_alt_stranded:True:4:False" \\
      --source "${hisat2_stringtie_default_unstranded}:hisat2_stringtie_default_unstranded:False:3:False" \\
      --source "${hisat2_stringtie_alt_unstranded}:hisat2_stringtie_alt_unstranded:False:2:False" \\
      --source "${long_reads_default}:long_reads_default:False:6:False" \\
      --source "${long_reads_alt}:long_reads_alt:False:5:False" \\
      --source "${flair_isoforms_gtf}:flair_isoforms:False:6:False" \\
      --source "${helixer_gff3}:helixer:False:4:False" \\
      -o transcript_inputs.tsv

    mikado configure \\
      --list transcript_inputs.tsv \\
      --reference "${genome}" \\
      --mode "${params.mikado_mode}" \\
      --scoring "${params.mikado_scoring}" \\
      mikado_configuration.yaml

    mikado prepare --json-conf mikado_configuration.yaml

    test -s mikado_prepared.fasta
    test -s mikado_prepared.gtf
    mikado --version 2>&1 | sed 's/^/  mikado: "/; s/\$/"/' | {
      printf '"%s":\\n' "${task.process}"
      cat
      printf '  container: "%s"\\n' "${task.container}"
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail
    seqid=\$(awk '/^>/ {print substr(\$1, 2); exit}' ${genome})
    printf "stub\\tmikado_stub\\tFalse\\t1\\tFalse\\n" > transcript_inputs.tsv
    printf "reference: %s\\n" "${genome}" > mikado_configuration.yaml
    printf ">mikado_stub_tx\\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\\n" > mikado_prepared.fasta
    printf "%s\\tMikado\\ttranscript\\t1\\t42\\t.\\t+\\t.\\tgene_id \\"mikado_stub_gene\\"; transcript_id \\"mikado_stub_tx\\";\\n" "\$seqid" > mikado_prepared.gtf
    printf '"%s":\n  mikado_prepare: "stub"\n' "${task.process}" > versions.yml
    """
}

process mikado_serialise {
  label 'process_transcriptome'
  tag "Mikado serialise with TransDecoder ORFs"
  container params.container_mikado

  input:
    path(config)
    path(prepared_fasta)
    path(transdecoder_bed)

  output:
    path "mikado.db", emit: database
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    if [[ "${params.run_mikado}" != "true" ]]; then
      : > mikado.db
      printf '"%s":\n  mikado_serialise: "skipped"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
      exit 0
    fi

    mikado serialise --json-conf "${config}" --orfs "${transdecoder_bed}" --procs ${task.cpus}
    test -s mikado.db
    mikado --version 2>&1 | sed 's/^/  mikado: "/; s/\$/"/' | {
      printf '"%s":\\n' "${task.process}"
      cat
      printf '  container: "%s"\\n' "${task.container}"
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "stub\\n" > mikado.db
    printf '"%s":\n  mikado_serialise: "stub"\n' "${task.process}" > versions.yml
    """
}

process mikado_pick {
  label 'process_transcriptome'
  tag "Mikado final annotation pick"
  container params.container_mikado
  publishDir "${params.output_dir}/final_annotations/mikado", mode: 'copy', saveAs: { filename ->
    if (filename in ['final_mikado_annotation.gff3', 'mikado.loci.gff3', 'mikado.subloci.gff3', 'versions.yml']) {
      return filename
    }
    return null
  }

  input:
    path(genome)
    path(config)
    path(database)

  output:
    path "final_mikado_annotation.gff3", emit: gff3
    path "mikado.loci.gff3", emit: loci_gff3
    path "mikado.subloci.gff3", emit: subloci_gff3
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    if [[ "${params.run_mikado}" != "true" ]]; then
      printf "##gff-version 3\\n" > final_mikado_annotation.gff3
      cp final_mikado_annotation.gff3 mikado.loci.gff3
      cp final_mikado_annotation.gff3 mikado.subloci.gff3
      printf '"%s":\n  mikado_pick: "skipped"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
      exit 0
    fi

    mikado pick --json-conf "${config}" --subloci-out mikado.subloci.gff3 --loci-out mikado.loci.gff3 --procs ${task.cpus}
    cp mikado.loci.gff3 final_mikado_annotation.gff3
    test -s final_mikado_annotation.gff3
    mikado --version 2>&1 | sed 's/^/  mikado: "/; s/\$/"/' | {
      printf '"%s":\\n' "${task.process}"
      cat
      printf '  container: "%s"\\n' "${task.container}"
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail
    seqid=\$(awk '/^>/ {print substr(\$1, 2); exit}' ${genome})
    end=\$(awk 'BEGIN{n=0} !/^>/ {n += length(\$0)} END{print n < 42 ? n : 42}' ${genome})
    cat > final_mikado_annotation.gff3 <<EOF
##gff-version 3
\${seqid}	Mikado	gene	1	\${end}	12.5	+	.	ID=mikado_stub_gene;mikado_score=12.5
\${seqid}	Mikado	mRNA	1	\${end}	12.5	+	.	ID=mikado_stub_tx;Parent=mikado_stub_gene;mikado_score=12.5
\${seqid}	Mikado	exon	1	\${end}	.	+	.	ID=mikado_stub_exon;Parent=mikado_stub_tx
\${seqid}	Mikado	CDS	1	\${end}	.	+	0	ID=mikado_stub_cds;Parent=mikado_stub_tx
EOF
    cp final_mikado_annotation.gff3 mikado.loci.gff3
    cp final_mikado_annotation.gff3 mikado.subloci.gff3
    printf '"%s":\n  mikado_pick: "stub"\n' "${task.process}" > versions.yml
    """
}

process final_annotation_sources_qc {
  label 'process_low'
  tag "AEGIS and Mikado final annotation source comparison"
  container params.container_python
  publishDir "${params.output_dir}/quality_report/final_annotation_sources", mode: 'copy'

  input:
    path(aegis_gff3)
    path(mikado_gff3)
    path(compare_script)

  output:
    path "final_annotation_sources.json", emit: json_report
    path "final_annotation_sources_mqc.tsv", emit: multiqc_tsv
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    python3 "${compare_script}" \\
      --aegis-gff3 "${aegis_gff3}" \\
      --mikado-gff3 "${mikado_gff3}" \\
      --json-report final_annotation_sources.json \\
      --multiqc-tsv final_annotation_sources_mqc.tsv
    printf '"%s":\n  compare_final_annotations: "python-stdlib"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    python3 "${compare_script}" \\
      --aegis-gff3 "${aegis_gff3}" \\
      --mikado-gff3 "${mikado_gff3}" \\
      --json-report final_annotation_sources.json \\
      --multiqc-tsv final_annotation_sources_mqc.tsv
    printf '"%s":\n  compare_final_annotations: "stub"\n' "${task.process}" > versions.yml
    """
}
