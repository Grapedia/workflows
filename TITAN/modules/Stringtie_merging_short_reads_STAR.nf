// Transcriptomes merging with StringTie
process Stringtie_merging_short_reads_STAR {

  tag "STAR/StringTie merging - short reads"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${params.output_dir}:/outputdir --volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work --volume ${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/STAR/stranded/:/StringTie_short_reads_STAR_stranded --volume ${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/STAR/unstranded/:/StringTie_short_reads_STAR_unstranded"
  cpus 4

  publishDir "${params.output_dir}/tmp", mode: 'copy'

  input:
    val(concat_star_stringtie_annot)

  output:
    path "merged_transcriptomes.STAR.short_reads.default_args.stranded.gtf", emit: default_args_stranded
    path "merged_transcriptomes.STAR.short_reads.alt_args.stranded.gtf", emit: alt_args_stranded
    path "merged_transcriptomes.STAR.short_reads.default_args.unstranded.gtf", optional: true, emit: default_args_unstranded
    path "merged_transcriptomes.STAR.short_reads.alt_args.unstranded.gtf", emit: alt_args_unstranded

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running StringTie merging on STAR/StringTie transcriptomes - separating stranded and unstranded samples."

    gtf_stranded_default=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads_STAR_stranded default)
    gtf_stranded_alt=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads_STAR_stranded alt)
    CMD="/scripts/Stringtie_merging.sh -o merged_transcriptomes.STAR.short_reads.default_args.stranded.gtf -g \${gtf_stranded_default}"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.STAR.short_reads.default_args.stranded.gtf -g \${gtf_stranded_default}
    CMD="/scripts/Stringtie_merging.sh -o merged_transcriptomes.STAR.short_reads.alt_args.stranded.gtf -g \${gtf_stranded_alt}"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.STAR.short_reads.alt_args.stranded.gtf -g \${gtf_stranded_alt}
    cp merged_transcriptomes.STAR.short_reads.default_args.stranded.gtf /outputdir/merged_star_stringtie_stranded_default.gtf
    cp merged_transcriptomes.STAR.short_reads.alt_args.stranded.gtf /outputdir/merged_star_stringtie_stranded_alt.gtf

    if [ -d "/StringTie_short_reads_STAR_unstranded" ] && [ "\$(ls -A /StringTie_short_reads_STAR_unstranded 2>/dev/null)" ]; then
      echo "[\$DATE] Running StringTie merging on STAR/StringTie transcriptomes - unstranded samples detected."
      gtf_unstranded_default=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads_STAR_unstranded default)
      gtf_unstranded_alt=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads_STAR_unstranded alt)
      CMD="/scripts/Stringtie_merging.sh -o merged_transcriptomes.STAR.short_reads.default_args.unstranded.gtf -g \${gtf_unstranded_default}"
      echo "[\$DATE] Executing: \$CMD"
      /scripts/Stringtie_merging.sh -o merged_transcriptomes.STAR.short_reads.default_args.unstranded.gtf -g \${gtf_unstranded_default}
      CMD="/scripts/Stringtie_merging.sh -o merged_transcriptomes.STAR.short_reads.alt_args.unstranded.gtf -g \${gtf_unstranded_alt}"
      echo "[\$DATE] Executing: \$CMD"
      /scripts/Stringtie_merging.sh -o merged_transcriptomes.STAR.short_reads.alt_args.unstranded.gtf -g \${gtf_unstranded_alt}
      cp merged_transcriptomes.STAR.short_reads.default_args.unstranded.gtf /outputdir/merged_star_stringtie_unstranded_default.gtf
      cp merged_transcriptomes.STAR.short_reads.alt_args.unstranded.gtf /outputdir/merged_star_stringtie_unstranded_alt.gtf
    fi
    """
}
