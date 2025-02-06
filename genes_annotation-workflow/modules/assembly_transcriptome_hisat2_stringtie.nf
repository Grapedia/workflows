// 4. Transcriptome assembly with StringTie
process assembly_transcriptome_hisat2_stringtie {

  tag "Hisat2/StringTie - short reads on ${sample_ID}"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work"
  publishDir { 
    if (strand_type == "unstranded") {
      return "$params.outdir/evidence_data/transcriptomes/StringTie/short_reads/HISAT2/unstranded"
    } else if (strand_type in ["stranded_forward", "stranded_reverse"]) {
      return "$params.outdir/evidence_data/transcriptomes/StringTie/short_reads/HISAT2/stranded"
    }
  }
  cpus 4

  input:
    tuple val(sample_ID), path(bam_file), val(strand_type)

  output:
    tuple val(sample_ID), path("*_transcriptome.gtf"), val(strand_type), emit: hisat2_stringtie_transcriptome_gtf
    tuple val(sample_ID), path("*_transcriptome.AltCommands.gtf"), val(strand_type), emit: hisat2_stringtie_alt_commands_gtf

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running transcriptome assembly with StringTie on HISAT2 $sample_ID"
    CMD="/scripts/Stringtie.sh -t ${task.cpus} -o ${sample_ID}_transcriptome.gtf -b ${bam_file} -r short"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/Stringtie.sh -t ${task.cpus} -o ${sample_ID}_transcriptome.gtf -b ${bam_file} -r short
    CMD="/scripts/Stringtie_AltCommands.sh -t ${task.cpus} -o ${sample_ID}_transcriptome.AltCommands.gtf -b ${bam_file} -r short"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/Stringtie_AltCommands.sh -t ${task.cpus} -o ${sample_ID}_transcriptome.AltCommands.gtf -b ${bam_file} -r short
    """
}
