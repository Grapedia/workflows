// Final step : gff3 merging with Aegis to create final GFF3 file
process aegis {

  tag "Generation of final GFF3 file with Aegis"
  // this image contains Aegis, Diamond and gffread
  container 'avelt/aegis:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work --volume $genome_path:/genome_path --volume ${projectDir}/data/protein_data:/protein_path --volume ${protein_samplesheet_path}:/protein_samplesheet_path"
  publishDir "$projectDir/FINAL_OUTPUT/aegis_final_GFF3_file/"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    val(protein_samplesheet_path)
    val(protein_samplesheet_filename)
    path(edta_masked_genome)
    path(augustus_gff)
    path(genemark_gtf)
    path(liftoff_annotations)
    path(stranded_default_args)
    path(stranded_alt_args)
    path(gffcompare_stranded)
    path(gffcompare_unstranded)
    path(unstranded_default_args, optional: true)
    path(unstranded_alt_args, optional: true)

  output:
    path("aegis_final_merged_annotations.gff3") emit: aegis_gff
    path("aegis_final_merged_annotations.pkl") emit: aegis_pkl
    path("*_unique_proteins.fasta") emit: aegis_proteins

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running Aegis to create final GFF3 file" >> ${params.logfile} 2>&1
    proteins=\$(/scripts/retrieve_proteins_for_braker.sh /protein_samplesheet_path/$protein_samplesheet_filename)
    # convert all GTF file in GFF3 format with gffread -E {gtf_file} -o- > {gff3_file}
    CMD_aegis_1="/scripts/Aegis1.py --genome_name New_assembly --genome_path /genome_path/$genome --augustus_path ${augustus_gff} --genemark_path ${genemark_gtf} --liftoff_path ${liftoff_annotations} --psiclass_stranded_STAR_path ${gffcompare_stranded} --stringtie_stranded_default_STAR_path ${stranded_default_args} --stringtie_stranded_AltCommands_STAR_path ${stranded_alt_args} --output_dir \$PWD --output_gff aegis_final_merged_annotations.gff3 --output_pickle aegis_final_merged_annotations.pkl
    if [ -s "${unstranded_default_args}" ]; then
        CMD_aegis_1="\$CMD_aegis_1 --stringtie_unstranded_default_STAR_path ${unstranded_default_args}"
    fi

    if [ -s "${unstranded_alt_args}" ]; then
        CMD_aegis_1="\$CMD_aegis_1 --stringtie_unstranded_AltCommands_STAR_path ${unstranded_alt_args}"
    fi
     echo "[\$DATE] Executing: \$CMD_aegis_1" >> ${params.logfile} 2>&1
     # /scripts/Aegis2.sh -q aegis_proteins.fasta -t ${task.cpus} -d \${proteins} -o \$PWD
     # /scripts/Aegis3.py --merged_annotation aegis_final_merged_annotations.pkl --hard_masked_genome ${edta_masked_genome} --blast_hits to_do --intermediate_annotation intermediate_annotation.pkl --final_annotation final_annotation.pkl --export_dir \$PWD --final_export_dir \$PWD --update --source_priority to_do
    """
}