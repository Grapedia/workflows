// 4. Transcriptome assembly with StringTie
process assembly_transcriptome_star_stringtie {

  tag "STAR/StringTie - short reads on ${sample_ID}"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work"
  publishDir { 
    if (strand_type == "unstranded") {
      return "${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/STAR/unstranded"
    } else if (strand_type in ["stranded_forward", "stranded_reverse"]) {
      return "${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/short_reads/STAR/stranded"
    }
  }

  cpus 4

  input:
    tuple val(sample_ID), path(bam_file), val(strand_type)

  output:
    tuple val(sample_ID), path("${sample_ID}_transcriptome.gtf"), path("${sample_ID}_transcriptome.AltCommands.gtf"), val(strand_type), emit: star_stringtie_transcriptomes

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running STAR/StringTie transcriptome assembly on $sample_ID"
    CMD="/scripts/Stringtie.sh -t ${task.cpus} -o ${sample_ID}_transcriptome.gtf -b ${bam_file} -r short"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/Stringtie.sh -t ${task.cpus} -o ${sample_ID}_transcriptome.gtf -b ${bam_file} -r short
    CMD="/scripts/Stringtie_AltCommands.sh -t ${task.cpus} -o ${sample_ID}_transcriptome.AltCommands.gtf -b ${bam_file} -r short"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/Stringtie_AltCommands.sh -t ${task.cpus} -o ${sample_ID}_transcriptome.AltCommands.gtf -b ${bam_file} -r short
    """

  stub:
    """
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_gene\\"; transcript_id \\"${sample_ID}_tx\\";\\n" > ${sample_ID}_transcriptome.gtf
    printf "chr1\\tStringTie\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"${sample_ID}_gene_alt\\"; transcript_id \\"${sample_ID}_tx_alt\\";\\n" > ${sample_ID}_transcriptome.AltCommands.gtf
    """
}
