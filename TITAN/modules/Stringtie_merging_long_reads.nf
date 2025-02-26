// Transcriptomes merging with StringTie
process Stringtie_merging_long_reads {

  tag "StringTie merging - long reads"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work --volume ${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/long_reads/:/StringTie_long_reads"
  cpus 4

  publishDir "${params.output_dir}/tmp", mode: 'copy'

  input:
    val(concat_minimap2_stringtie_annot)

  output:
    path "merged_transcriptomes.minimap2.long_reads.default_args.gtf", emit: default_args_gff
    path "merged_transcriptomes.minimap2.long_reads.alt_args.gtf", emit: alt_args_gff

  when:
  params.use_long_reads

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running Minimap2/StringTie merging - long reads transcriptome assemblies"
    gtf_default=\$(/scripts/retrieve_path_gtf.sh /StringTie_long_reads default)
    gtf_alt=\$(/scripts/retrieve_path_gtf.sh /StringTie_long_reads alt)
    CMD="/scripts/Stringtie_merging.sh -o merged_transcriptomes.minimap2.long_reads.default_args.gtf -g \${gtf_default}"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.minimap2.long_reads.default_args.gtf -g \${gtf_default}
    CMD="/scripts/Stringtie_merging.sh -o merged_transcriptomes.minimap2.long_reads.alt_args.gtf -g \${gtf_alt}"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.minimap2.long_reads.alt_args.gtf -g \${gtf_alt}
    cp merged_transcriptomes.minimap2.long_reads.default_args.gtf ${params.output_dir}/merged_minimap2_stringtie_long_reads_default.gtf
    cp merged_transcriptomes.minimap2.long_reads.alt_args.gtf ${params.output_dir}/merged_minimap2_stringtie_long_reads_alt.gtf
    """
}
