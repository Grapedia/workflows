// We check that if the samples are of type “FASTQ” they exist in the folder data/RNAseq_data.
// If the samples are SRA, we'll download the .sra file, convert it to fastq and gzip the fastq file.
// Once done, the RNAseq samples are ready for transcriptome assembly.

process prepare_RNAseq_fastq_files_short {
  tag "prepare_RNAseq_fastq_files on $sample_ID"
  container 'quay.io/biocontainers/sra-tools:3.1.1--h4304569_0'
  containerOptions "--volume ${projectDir}/data/RNAseq_data:/RNAseq_data"
  echo true
  // Sometimes the SRA download encounters an error, so we retry the process until the download works.
  errorStrategy 'retry'

  input:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout)

  output:
  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout)

  script:
  """
  DATE=\$(date "+%Y-%m-%d %H:%M:%S")
  echo "[\$DATE] Running prepare_RNAseq_fastq_files_short on $sample_ID"

  if [[ $SRA_or_FASTQ == "SRA" ]]
  then
    if [[ $library_layout == "paired" ]]
    then
      if ls /RNAseq_data/${sample_ID}*.fastq.gz 1> /dev/null 2>&1
      then
        echo "[\$DATE] Sample $sample_ID already processed. Skipping."
        continue
      else
        echo "[\$DATE] Downloading SRA sample: $sample_ID - paired-end sample"
        prefetch --force all -O /RNAseq_data/ $sample_ID
        echo "[\$DATE] Converting .sra to FASTQ for $sample_ID"
        fastq-dump --outdir /RNAseq_data/ --split-files $sample_ID
        rm -R /RNAseq_data/$sample_ID
        echo "[\$DATE] Compressing FASTQ files for $sample_ID"
        gzip /RNAseq_data/${sample_ID}_1.fastq
        gzip /RNAseq_data/${sample_ID}_2.fastq
      fi
    elif [[ $library_layout == "single" ]]
    then
      if ls /RNAseq_data/${sample_ID}.fastq.gz 1> /dev/null 2>&1
      then
        echo "[\$DATE] Sample $sample_ID already processed. Skipping."
        continue
      else
        echo "[\$DATE] Downloading SRA sample: $sample_ID - single-end sample"
        prefetch --force all -O /RNAseq_data/ $sample_ID
        echo "[\$DATE] Converting .sra to FASTQ for $sample_ID"
        fastq-dump --outdir /RNAseq_data/ $sample_ID
        rm -R /RNAseq_data/$sample_ID
        echo "[\$DATE] Compressing FASTQ file for $sample_ID"
        gzip /RNAseq_data/${sample_ID}.fastq
      fi
    else
      echo "[\$DATE] WARNING : $library_layout is not equal to paired or single !"
    fi
  elif [[ $SRA_or_FASTQ == "FASTQ" ]]
  then
    if ls /RNAseq_data/$sample_ID*fastq.gz 1> /dev/null 2>&1
    then
      echo "[\$DATE] FASTQ sample $sample_ID already exists. Skipping."
      continue
    else
      echo "[\$DATE] WARNING: ${projectDir}/data/RNAseq_data/$sample_ID*fastq.gz does not exist!"
    fi
  else
    echo "[\$DATE] ERROR: \$SRA_or_FASTQ is not equal to SRA or FASTQ!"
  fi
  """
}
