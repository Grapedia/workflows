process pasa {

  tag "Executing PASA UTR annotation with {threads} threads on the following GFF3: ${gff3_file}"
  container 'quay.io/biocontainers/pasa:2.4.1--he1b5a44_0'
  containerOptions "--volume ${genome_path}:/genome_path --volume ${pasa_config_path}:/pasa_config_path"
  publishDir "$projectDir/FINAL_OUTPUT"
  cpus 4

  input:
    val(gff3_file)
    val(genome_path)
    val(genome)
    val(transcriptome)
    val(pasa_config_path)
    val(pasa_config_filename)

  output:
    file("annotations.filtered.UTRs.gff3")
    file("pasa.gene_structures_post_PASA_updates.*.gff3")

  script:
    """
    /usr/local/opt/pasa-2.4.1/Launch_PASA_pipeline.pl -c /pasa_config_path/${pasa_config_filename} --annot_compare -r -C -R -g /genome_path/${genome} -L --annot_compare --annots ${gff3_file} \
    -t ${transcriptome} \
    --stringent_alignment_overlap 30 --TRANSDECODER --MAX_INTRON_LENGTH 3000 --CPU ${task.cpus} --ALIGNERS blat
    """
}
