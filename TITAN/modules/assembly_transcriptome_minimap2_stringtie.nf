// 4. Transcriptome assembly with StringTie
process assembly_transcriptome_minimap2_stringtie {

  tag "Minimap2/StringTie - long reads"
  container 'avelt/stringtie:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work"
  publishDir "${params.output_dir}/intermediate_files/evidence_data/transcriptomes/StringTie/long_reads"
  cpus 4

  input:
    tuple val(sample_ID), path(bam_file)

  output:
    tuple val(sample_ID), path("*_transcriptome.gtf"), emit: minimap2_stringtie_transcriptome_gtf
    tuple val(sample_ID), path("*_transcriptome.AltCommands.gtf"), emit: minimap2_stringtie_alt_commands_gtf

  when:
  params.use_long_reads

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running Minimap2/StringTie - long reads transcriptome assembly on $sample_ID"
    CMD="/scripts/Stringtie.sh -t ${task.cpus} -o ${sample_ID}_transcriptome.gtf -b ${bam_file} -r long"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/Stringtie.sh -t ${task.cpus} -o ${sample_ID}_transcriptome.gtf -b ${bam_file} -r long
    CMD="/scripts/Stringtie_AltCommands.sh -o ${sample_ID}_transcriptome.AltCommands.gtf -b ${bam_file} -r long"
    echo "[\$DATE] Executing: \$CMD"
    /scripts/Stringtie_AltCommands.sh -o ${sample_ID}_transcriptome.AltCommands.gtf -b ${bam_file} -r long
    """
}
