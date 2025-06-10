// Final step : gff3 merging with Aegis to create final GFF3 file
process aegis_long_reads {

  tag "Generation of final GFF3 file with Aegis"
  // this image contains Aegis, Diamond and gffread
  container 'avelt/aegis:v2025_05_20'
  containerOptions "--volume ${projectDir}:/projectDir --volume ${projectDir}/scripts/:/scripts --volume ${projectDir}/work:/work --volume $genome_path:/genome_path --volume ${projectDir}/data/protein_data:/protein_path --volume ${protein_samplesheet_path}:/protein_samplesheet_path"
  publishDir "${params.output_dir}/aegis_outputs", mode: 'copy'
  cpus 1

  input:
    val(genome_path)
    val(genome)
    val(protein_samplesheet_path)
    val(protein_samplesheet_filename)
    path(edta_masked_genome)
    path(augustus_gff)
    path(genemark_gtf)
    path(liftoff_annotations)
    path(long_reads_default_args)
    path(long_reads_alt_args)
    path(stranded_default_args)
    path(stranded_alt_args)
    path(gffcompare_stranded)
    path(gffcompare_unstranded)
    path(unstranded_default_args)
    path(unstranded_alt_args)

  output:
    path "final_annotation.gff3", emit: aegis_gff
    path "final_annotation_proteins_all.fasta", emit: aegis_proteins_all
    path "final_annotation_proteins_main.fasta", emit: aegis_proteins_main

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running Aegis to create final GFF3 file"

    mkdir -p /projectDir/tmp
    export TMPDIR=/projectDir/tmp

    mapfile -t proteins < <(/scripts/retrieve_proteins_for_aegis.sh /protein_samplesheet_path/$protein_samplesheet_filename)

    # create a bash array with the name of the samplesheet organism as the first value and the path to the corresponding fasta file as the second value
    declare -A PROTEIN_MAP
    declare -a PROTEIN_KEYS

    for line in "\${proteins[@]}"; do
        clean_line=\$(echo "\$line" | tr -d '\r' | awk '{\$1=\$1; print}')

        IFS=' ' read -r name path <<< "\$clean_line"

        PROTEIN_MAP["\$name"]="\$path"
        PROTEIN_KEYS+=("\$name")
    done

    # here we create the protein_paths parameters for Aegis2.sh of type : /path/to/arabidopsis.proteins.fasta,/path/to/viridiplantae.proteins.fasta
    protein_paths=""
    for key in "\${!PROTEIN_MAP[@]}"; do
        if [[ -z "\$protein_paths" ]]; then
            protein_paths="\${PROTEIN_MAP[\$key]}"
        else
            protein_paths+=",\${PROTEIN_MAP[\$key]}"
        fi
    done

    # here we create the source_priority parameters for Aegis3.py of type : ['Araport', 'Viridiplantae', 'Eudicots']
    protein_names="["
    first_entry=true
    for key in "\${PROTEIN_KEYS[@]}"; do
        if \$first_entry; then
            protein_names+="'\$key'"
            first_entry=false
        else
            protein_names+=", '\$key'"
        fi
    done
    protein_names+="]"

    # gtf to gff3 conversion
    gffread -E ${genemark_gtf} -o- > BRAKER3_genemark.gff3
    gffread -E ${stranded_default_args} -o- > stringtie_stranded_default_STAR.gff3
    gffread -E ${stranded_alt_args} -o- > stringtie_stranded_AltCommands_STAR.gff3
    gffread -E ${gffcompare_stranded} -o- > psiclass_stranded_STAR.gff3
    gffread -E ${long_reads_default_args} -o- > stringtie_Isoseq_default.gff3
    gffread -E ${long_reads_alt_args} -o- > stringtie_Isoseq_AltCommands.gff3

    if [ "\$(basename ${gffcompare_unstranded})" != "dev_null3" ]; then
        gffread -E ${gffcompare_unstranded} -o- > psiclass_unstranded_STAR.gff3
    fi
    if [ "\$(basename ${unstranded_default_args})" != "dev_null1" ]; then
        gffread -E ${unstranded_default_args} -o- > stringtie_unstranded_default_STAR.gff3
    fi
    if [ "\$(basename ${unstranded_alt_args})" != "dev_null2" ]; then
        gffread -E ${unstranded_alt_args} -o- > stringtie_unstranded_AltCommands_STAR.gff3
    fi

    # Extract chromosomes names from fasta file
    chromosomes=\$(grep '^>' /genome_path/$genome | cut -d' ' -f1 | sed 's/>//' | uniq)

    for chrom in \$chromosomes
    do
        # here we create the path to Diamond results, of type :
        # Viridiplantae=viridiplantae_vs_proteins_assembly.diamond Eudicots=eudicots_vs_proteins_assembly.diamond
        diamond_paths=""
        for key in "\${!PROTEIN_MAP[@]}"; do
            file_name=\$(basename "\${PROTEIN_MAP[\$key]}")
            file_path="\${chrom}/\${file_name}_vs_assembly.diamond"
            if [[ -z "\$diamond_paths" ]]; then
                diamond_paths="\$key=\$file_path"
            else
                diamond_paths+=" \$key=\$file_path"
            fi
        done

        echo "[\$(date '+%Y-%m-%d %H:%M:%S')] Aegis on chromosome : \$chrom"

        CMD_aegis_1="/scripts/Aegis1.py --genome_name New_assembly --genome_path /genome_path/$genome --augustus_path ${augustus_gff} --genemark_path BRAKER3_genemark.gff3 --liftoff_path ${liftoff_annotations} --psiclass_stranded_STAR_path psiclass_stranded_STAR.gff3 --stringtie_stranded_default_STAR_path stringtie_stranded_default_STAR.gff3 --stringtie_stranded_AltCommands_STAR_path stringtie_stranded_AltCommands_STAR.gff3 --stringtie_Isoseq_default_path stringtie_Isoseq_default.gff3 --stringtie_Isoseq_AltCommands_path stringtie_Isoseq_AltCommands.gff3 --output_dir \$PWD/\${chrom} --output_gff aegis_final_merged_annotations_\${chrom}.gff3 --output_pickle aegis_final_merged_annotations_\${chrom}.pkl --chromosome \$chrom"
        if [ "\$(basename ${gffcompare_unstranded})" != "dev_null3" ]; then
            CMD_aegis_1="\$CMD_aegis_1 --psiclass_unstranded_STAR_path psiclass_unstranded_STAR.gff3"
        fi
        if [ "\$(basename ${unstranded_default_args})" != "dev_null1" ]; then
            CMD_aegis_1="\$CMD_aegis_1 --stringtie_unstranded_default_STAR_path stringtie_unstranded_default_STAR.gff3"
        fi
        if [ "\$(basename ${unstranded_alt_args})" != "dev_null2" ]; then
            CMD_aegis_1="\$CMD_aegis_1 --stringtie_unstranded_AltCommands_STAR_path stringtie_unstranded_AltCommands_STAR.gff3"
        fi
         echo "[\$DATE] Executing: \$CMD_aegis_1"
         eval "\$CMD_aegis_1"
         echo "[\$DATE] Executing: /scripts/Aegis2.sh -q \$PWD/\${chrom}/merged_annotation_\${chrom}_unique_proteins.fasta -t ${task.cpus} -d \$protein_paths -o \$PWD/\${chrom}"
         /scripts/Aegis2.sh -q \$PWD/\${chrom}/merged_annotation_\${chrom}_unique_proteins.fasta -t ${task.cpus} -d \$protein_paths -o \$PWD/\${chrom}
         echo "[\$DATE] Executing: /scripts/Aegis3.py --merged_annotation \$PWD/\${chrom}/aegis_final_merged_annotations_\${chrom}.pkl --hard_masked_genome ${edta_masked_genome} --diamond_hits \${diamond_paths} --intermediate_annotation intermediate_annotation_\${chrom}.pkl --final_annotation final_annotation_\${chrom}.pkl --export_dir \$PWD/\${chrom} --final_export_dir \$PWD/\${chrom} --update --chromosome \$chrom --source_priority \$protein_names"
         /scripts/Aegis3.py --merged_annotation \$PWD/\${chrom}/aegis_final_merged_annotations_\${chrom}.pkl --hard_masked_genome ${edta_masked_genome} --diamond_hits \${diamond_paths} --intermediate_annotation intermediate_annotation_\${chrom}.pkl --final_annotation final_annotation_\${chrom}.pkl --export_dir \$PWD/\${chrom} --final_export_dir \$PWD/\${chrom} --update --chromosome \$chrom --source_priority \$protein_names
    done

    # final concatenate of all final GFF3 files
    gff3_final="\$PWD/final_annotation.gff3"

    echo "##gff-version 3" > \$gff3_final

    for gff3_chrom in `ls \$PWD/*/final_annotation_*.gff3`
    do
        tail -n +4 \$gff3_chrom >> \$gff3_final
    done

    # final concatenate of all final protein files
    proteins_final_all="\$PWD/final_annotation_proteins_all.fasta"
    proteins_final_main="\$PWD/final_annotation_proteins_main.fasta"

    > "\$proteins_final_all"

    for proteins_chrom_all in `ls \$PWD/*/final_annotation_*_proteins_p_id_all.fasta`
    do
        cat \$proteins_chrom_all >> \$proteins_final_all
    done

    > "\$proteins_final_main"

    for proteins_chrom_main in `ls \$PWD/*/final_annotation_*_proteins_p_id_main.fasta`
    do
        cat \$proteins_chrom_main >> \$proteins_final_main
    done
    """
}
