// Transcriptomes merging with StringTie
process Stringtie_merging_long_reads {
  label 'process_merge'

  tag "StringTie merging - long reads"
  container params.container_stringtie
  publishDir "${params.output_dir}/tmp", mode: 'copy'
  publishDir "${params.output_dir}", mode: 'copy', saveAs: { filename ->
    if (filename == 'merged_transcriptomes.minimap2.long_reads.default_args.gtf') {
      return 'merged_minimap2_stringtie_long_reads_default.gtf'
    }
    if (filename == 'merged_transcriptomes.minimap2.long_reads.alt_args.gtf') {
      return 'merged_minimap2_stringtie_long_reads_alt.gtf'
    }
    return null
  }

  input:
    path(default_gtfs)
    path(alt_gtfs)

  output:
    path "merged_transcriptomes.minimap2.long_reads.default_args.gtf", emit: default_args_gff
    path "merged_transcriptomes.minimap2.long_reads.alt_args.gtf", emit: alt_args_gff

  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running Minimap2/StringTie merging - long reads transcriptome assemblies"
    printf '%s\\n' ${default_gtfs} > default_gtfs.txt
    printf '%s\\n' ${alt_gtfs} > alt_gtfs.txt
    stringtie --merge -o merged_transcriptomes.minimap2.long_reads.default_args.gtf default_gtfs.txt
    stringtie --merge -o merged_transcriptomes.minimap2.long_reads.alt_args.gtf alt_gtfs.txt
    """

  stub:
    """
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"long_merged_gene\\"; transcript_id \\"long_merged_tx\\";\\n" > merged_transcriptomes.minimap2.long_reads.default_args.gtf
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"long_merged_alt_gene\\"; transcript_id \\"long_merged_alt_tx\\";\\n" > merged_transcriptomes.minimap2.long_reads.alt_args.gtf
    """
}
