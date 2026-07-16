// We check that if the samples are of type FASTQ they exist in params.RNAseq_data_dir.
// If the samples are SRA, we'll download the .sra file, convert it to fastq and gzip the fastq file.
// Once done, the RNAseq samples are ready for transcriptome assembly.

process prepare_RNAseq_fastq_files_long {
  tag "prepare_RNAseq_fastq_files on $sample_ID"
  container params.container_sra_tools
  debug true
  // Sometimes the SRA download encounters an error, so we retry the process until the download works.
  errorStrategy 'retry'

  input:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path(local_reads)

  output:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path("${sample_ID}*", includeInputs: true), emit : prepared_fastqs

  script:
  """
  set -euo pipefail
  DATE=\$(date "+%Y-%m-%d %H:%M:%S")
  echo "[\$DATE] Preparing long RNA-seq files for $sample_ID"

  if [[ $SRA_or_FASTQ == "SRA" ]]
  then
    prefetch --force all --max-size 100G -O . $sample_ID
    fastq-dump --outdir . $sample_ID
    rm -rf "$sample_ID"
    gzip ${sample_ID}*.fastq
  elif [[ $SRA_or_FASTQ == "FASTQ" ]]
  then
    test -s ${sample_ID}.fastq.gz
  elif [[ $SRA_or_FASTQ == "FASTA" ]]
  then
    test -s ${sample_ID}.fasta
  else
    echo "[\$DATE] ERROR: \$SRA_or_FASTQ is not equal to SRA, FASTQ or FASTA" >&2
    exit 1
  fi
  """

  stub:
  """
  rm -f ${sample_ID}*
  printf "@stub\\nACGT\\n+\\n!!!!\\n" | gzip -c > ${sample_ID}.fastq.gz
  """
}
