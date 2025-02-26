// Transcriptomes merging with StringTie
process Stringtie_merging_short_reads_hisat2 {

  tag "Hisat2/StringTie merging - short reads"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work --volume ${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/HISAT2/stranded/:/StringTie_short_reads_HISAT2_stranded --volume ${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/HISAT2/unstranded/:/StringTie_short_reads_HISAT2_unstranded"
  publishDir "${params.output_dir}", mode: 'copy'
  cpus 4

  input:
    val(concat_star_stringtie_annot)

  output:
    file("*.gtf")

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running StringTie merging on HISAT2/StringTie transcriptomes - separating stranded and unstranded samples."

    gtf_stranded_default=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads_HISAT2_stranded default)
    gtf_stranded_alt=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads_HISAT2_stranded alt)
    CMD="/scripts/Stringtie_merging.sh -o merged_transcriptomes.hisat2.short_reads.default_args.stranded.gtf -g \${gtf_stranded_default}"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.hisat2.short_reads.default_args.stranded.gtf -g \${gtf_stranded_default}
    CMD="/scripts/Stringtie_merging.sh -o merged_transcriptomes.hisat2.short_reads.alt_args.stranded.gtf -g \${gtf_stranded_alt}"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/Stringtie_merging.sh -o merged_transcriptomes.hisat2.short_reads.alt_args.stranded.gtf -g \${gtf_stranded_alt}

    if [ -d "/StringTie_short_reads_HISAT2_unstranded" ] && [ "\$(ls -A /StringTie_short_reads_HISAT2_unstranded 2>/dev/null)" ]; then
        echo "[\$DATE] Running StringTie merging on HISAT2/StringTie transcriptomes - unstranded samples detected."
        gtf_unstranded_default=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads_HISAT2_unstranded default)
        gtf_unstranded_alt=\$(/scripts/retrieve_path_gtf.sh /StringTie_short_reads_HISAT2_unstranded alt)
        CMD="/scripts/Stringtie_merging.sh -o merged_transcriptomes.hisat2.short_reads.default_args.unstranded.gtf -g \${gtf_unstranded_default}"
        echo "[\$DATE] Executing: \$CMD"
        /scripts/Stringtie_merging.sh -o merged_transcriptomes.hisat2.short_reads.default_args.unstranded.gtf -g \${gtf_unstranded_default}
        CMD="/scripts/Stringtie_merging.sh -o merged_transcriptomes.hisat2.short_reads.alt_args.unstranded.gtf -g \${gtf_unstranded_alt}"
        echo "[\$DATE] Executing: \$CMD"
        /scripts/Stringtie_merging.sh -o merged_transcriptomes.hisat2.short_reads.alt_args.unstranded.gtf -g \${gtf_unstranded_alt}
    fi
    """
}
