// We check that if the samples are of type FASTQ they exist in params.RNAseq_data_dir.
// If the samples are SRA, we'll download the .sra file, convert it to fastq and gzip the fastq file.
// Once done, the RNAseq samples are ready for transcriptome assembly.

process prepare_RNAseq_fastq_files_short {
  label 'process_low'
  tag "prepare_RNAseq_fastq_files on $sample_ID"
  container params.container_sra_tools
  debug true
  // Sometimes the SRA download encounters an error, so we retry the process until the download works.
  errorStrategy 'retry'

  input:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path(local_reads)

  output:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path("${sample_ID}*.fastq.gz", includeInputs: true), emit: prepared_fastqs

  script:
  """
  set -euo pipefail
  DATE=\$(date "+%Y-%m-%d %H:%M:%S")
  echo "[\$DATE] Preparing short RNA-seq files for $sample_ID"

  if [[ $SRA_or_FASTQ == "SRA" ]]
  then
    if [[ $library_layout == "paired" ]]
    then
      echo "[\$DATE] Downloading SRA paired-end sample: $sample_ID"
      prefetch --force all --max-size 100G -O . $sample_ID
      fastq-dump --outdir . --split-files $sample_ID
      rm -rf "$sample_ID"
      gzip ${sample_ID}_1.fastq ${sample_ID}_2.fastq
    elif [[ $library_layout == "single" ]]
    then
      echo "[\$DATE] Downloading SRA single-end sample: $sample_ID"
      prefetch --force all --max-size 100G -O . $sample_ID
      fastq-dump --outdir . $sample_ID
      rm -rf "$sample_ID"
      gzip ${sample_ID}.fastq
    else
      echo "[\$DATE] ERROR: $library_layout is not equal to paired or single" >&2
      exit 1
    fi
  elif [[ $SRA_or_FASTQ == "FASTQ" ]]
  then
    if [[ $library_layout == "paired" ]]; then
      test -s ${sample_ID}_1.fastq.gz
      test -s ${sample_ID}_2.fastq.gz
    elif [[ $library_layout == "single" ]]; then
      test -s ${sample_ID}.fastq.gz
    else
      echo "[\$DATE] ERROR: $library_layout is not equal to paired or single" >&2
      exit 1
    fi
  else
    echo "[\$DATE] ERROR: \$SRA_or_FASTQ is not equal to SRA or FASTQ" >&2
    exit 1
  fi
  """

  stub:
  """
  rm -f ${sample_ID}*.fastq.gz
  if [[ $library_layout == "paired" ]]
  then
    printf "@stub/1\\nACGT\\n+\\n!!!!\\n" | gzip -c > ${sample_ID}_1.fastq.gz
    printf "@stub/2\\nTGCA\\n+\\n!!!!\\n" | gzip -c > ${sample_ID}_2.fastq.gz
  else
    printf "@stub\\nACGT\\n+\\n!!!!\\n" | gzip -c > ${sample_ID}.fastq.gz
  fi
  """
}
