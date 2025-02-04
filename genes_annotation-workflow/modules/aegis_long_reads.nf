// Final step : gff3 merging with Aegis to create final GFF3 file
process aegis_long_reads {

  tag "Generation of final GFF3 file with Aegis"
  // this image contains Aegis, Diamond and gffread
  container 'avelt/aegis:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work --volume $genome_path:/genome_path"
  publishDir "$projectDir/FINAL_OUTPUT/aegis_final_GFF3_file/"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    path augustus_gff
    path genemark_gtf
    path liftoff_annotations
    path long_reads_default_args
    path long_reads_alt_args
    path stranded_default_args
    path stranded_alt_args
    path gffcompare_stranded
    path gffcompare_unstranded
    path unstranded_default_args, optional: true
    path unstranded_alt_args, optional: true

  output:
    path("aegis_final_merged_annotations.gff3") emit: aegis_gff
    path("aegis_final_merged_annotations.pkl") emit: aegis_pkl
    path("*_unique_proteins.fasta") emit: aegis_proteins

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running Aegis to create final GFF3 file" >> ${params.logfile} 2>&1
    # convert all GTF file in GFF3 format with gffread -E {gtf_file} -o- > {gff3_file}
    CMD_aegis_1="/scripts/Aegis1.py --genome_name New_assembly --genome_path /genome_path/$genome --augustus_path ${augustus_gff} --genemark_path ${genemark_gtf} --liftoff_path ${liftoff_annotations} --psiclass_stranded_STAR_path ${gffcompare_stranded} --stringtie_stranded_default_STAR_path ${stranded_default_args} --stringtie_stranded_AltCommands_STAR_path ${stranded_alt_args} --stringtie_Isoseq_default_path ${long_reads_default_args} --stringtie_Isoseq_AltCommands_path ${long_reads_alt_args} --output_dir \$PWD --output_gff aegis_final_merged_annotations.gff3 --output_pickle aegis_final_merged_annotations.pkl
    if [ -s "${unstranded_default_args}" ]; then
        CMD_aegis_1="\$CMD_aegis_1 --stringtie_unstranded_default_STAR_path ${unstranded_default_args}"
    fi

    if [ -s "${unstranded_alt_args}" ]; then
        CMD_aegis_1="\$CMD_aegis_1 --stringtie_unstranded_AltCommands_STAR_path ${unstranded_alt_args}"
    fi
     echo "[\$DATE] Executing: \$CMD_aegis_1" >> ${params.logfile} 2>&1
     # /scripts/Aegis2.sh
     # /scripts/Aegis3.py
    """
}
