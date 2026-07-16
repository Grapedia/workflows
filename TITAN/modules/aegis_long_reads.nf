process aegis_long_reads {
  label 'process_aegis'
  tag "AEGIS merge of short- and long-read evidence"

  container "${params.aegis_container}"
  publishDir "${params.output_dir}/aegis_outputs", mode: 'copy'
  input:
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

  output:
    path "final_annotation.gff3", emit: aegis_gff
    path "final_annotation_proteins_all.fasta", emit: aegis_proteins_all
    path "final_annotation_proteins_main.fasta", emit: aegis_proteins_main
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running AEGIS merge on short- and long-read evidence"

    merge_inputs=(
      "${liftoff_annotations}"
      "${augustus_gff}"
      "${genemark_gtf}"
      "${egapx_gff3}"
      "${stranded_default_args}"
      "${stranded_alt_args}"
      "${gffcompare_stranded}"
      "${long_reads_default_args}"
      "${long_reads_alt_args}"
    )

    if [[ -s "${gffcompare_unstranded}" && "\$(basename "${gffcompare_unstranded}")" != dev_null* ]]; then
      merge_inputs+=("${gffcompare_unstranded}")
    fi
    if [[ -s "${unstranded_default_args}" && "\$(basename "${unstranded_default_args}")" != dev_null* ]]; then
      merge_inputs+=("${unstranded_default_args}")
    fi
    if [[ -s "${unstranded_alt_args}" && "\$(basename "${unstranded_alt_args}")" != dev_null* ]]; then
      merge_inputs+=("${unstranded_alt_args}")
    fi

    printf "[\$DATE] AEGIS merge inputs:\\n"
    printf '  %s\\n' "\${merge_inputs[@]}"

    aegis merge -d aegis_merge -o final_annotation "\${merge_inputs[@]}"
    cp aegis_merge/final_annotation.gff3 final_annotation.gff3

    aegis extract -f protein -m all -d aegis_proteins_all final_annotation.gff3 "${edta_masked_genome}"
    aegis extract -f protein -m main -d aegis_proteins_main final_annotation.gff3 "${edta_masked_genome}"

    all_proteins=\$(find aegis_proteins_all -type f -name '*proteins*all*.fasta' -print -quit 2>/dev/null || true)
    main_proteins=\$(find aegis_proteins_main -type f -name '*proteins*main*.fasta' -print -quit 2>/dev/null || true)

    if [[ -n "\$all_proteins" ]]; then
      cp "\$all_proteins" final_annotation_proteins_all.fasta
    else
      echo "AEGIS did not produce an all-proteins FASTA" >&2
      exit 1
    fi

    if [[ -n "\$main_proteins" ]]; then
      cp "\$main_proteins" final_annotation_proteins_main.fasta
    else
      echo "AEGIS did not produce a main-proteins FASTA" >&2
      exit 1
    fi

    printf '"%s":\\n  aegis: "%s"\\n  aegis_container: "%s"\\n' \\
      "${task.process}" "${params.aegis_version}" "${params.aegis_container}" > versions.yml
    """

  stub:
    """
    printf "##gff-version 3\\nchr1\\tAegis\\tgene\\t1\\t10\\t.\\t+\\t.\\tID=aegis_stub_gene\\n" > final_annotation.gff3
    printf ">aegis_stub_all\\nM\\n" > final_annotation_proteins_all.fasta
    printf ">aegis_stub_main\\nM\\n" > final_annotation_proteins_main.fasta
    printf '"%s":\\n  aegis: "%s"\\n  aegis_container: "%s"\\n' \\
      "${task.process}" "${params.aegis_version}" "${params.aegis_container}" > versions.yml
    """
}
