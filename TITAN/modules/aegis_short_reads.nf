process aegis_short_reads {
    tag "Generation of final GFF3 file with Aegis"

    container 'avelt/aegis:latest'
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
        path(stranded_default_args)
        path(stranded_alt_args)
        path(gffcompare_stranded)
        path(gffcompare_unstranded)
        path(unstranded_default_args)
        path(unstranded_alt_args)

    output:
        path "final_annotation.gff3", emit: aegis_gff
        path "final_annotation.pkl", emit: aegis_pkl
        path "final_annotation_main_proteins.fasta", emit: aegis_proteins_main
        path "final_annotation_all_proteins.fasta", emit: aegis_proteins_all

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

        # here we create the path to Diamond results, of type : 
        # Viridiplantae=/path/to/viridiplantae_vs_proteins_assembly.diamond Eudicots=/path/to/eudicots_vs_proteins_assembly.diamond
        diamond_paths=""
        for key in "\${!PROTEIN_MAP[@]}"; do
            file_name=\$(basename "\${PROTEIN_MAP[\$key]}")
            file_path="\${file_name}_vs_assembly.diamond"

            if [[ -z "\$diamond_paths" ]]; then
                diamond_paths="\$key=\$file_path"
            else
                diamond_paths+=" \$key=\$file_path"
            fi
        done

        # gtf to gff3 conversion
        gffread -E ${genemark_gtf} -o- > BRAKER3_genemark.gff3
        gffread -E ${stranded_default_args} -o- > stringtie_stranded_default_STAR.gff3
        gffread -E ${stranded_alt_args} -o- > stringtie_stranded_AltCommands_STAR.gff3
        gffread -E ${gffcompare_stranded} -o- > psiclass_stranded_STAR.gff3
        if [ -s "${gffcompare_unstranded}" ]; then
            gffread -E ${gffcompare_unstranded} -o- > psiclass_unstranded_STAR.gff3
        fi
        if [ -s "${unstranded_default_args}" ]; then
            gffread -E ${unstranded_default_args} -o- > stringtie_unstranded_default_STAR.gff3
        fi
        if [ -s "${unstranded_alt_args}" ]; then
            gffread -E ${unstranded_alt_args} -o- > stringtie_unstranded_AltCommands_STAR.gff3
        fi

        CMD_aegis_1="/scripts/Aegis1.py --genome_name New_assembly --genome_path /genome_path/$genome --augustus_path ${augustus_gff} --genemark_path BRAKER3_genemark.gff3 --liftoff_path ${liftoff_annotations} --psiclass_stranded_STAR_path psiclass_stranded_STAR.gff3 --stringtie_stranded_default_STAR_path stringtie_stranded_default_STAR.gff3 --stringtie_stranded_AltCommands_STAR_path stringtie_stranded_AltCommands_STAR.gff3 --output_dir \$PWD --output_gff aegis_final_merged_annotations.gff3 --output_pickle aegis_final_merged_annotations.pkl"
        if [ -s "${gffcompare_unstranded}" ]; then
            CMD_aegis_1="\$CMD_aegis_1 --psiclass_unstranded_STAR_path psiclass_unstranded_STAR.gff3"
        fi
        if [ -s "${unstranded_default_args}" ]; then
            CMD_aegis_1="\$CMD_aegis_1 --stringtie_unstranded_default_STAR_path stringtie_unstranded_default_STAR.gff3"
        fi
        if [ -s "${unstranded_alt_args}" ]; then
            CMD_aegis_1="\$CMD_aegis_1 --stringtie_unstranded_AltCommands_STAR_path stringtie_unstranded_AltCommands_STAR.gff3"
        fi

         echo "[\$DATE] Executing: \$CMD_aegis_1"
         eval "\$CMD_aegis_1"
         echo "[\$DATE] Executing: /scripts/Aegis2.sh -q merge_annotation_unique_proteins.fasta -t ${task.cpus} -d \$protein_paths -o \$PWD"
         /scripts/Aegis2.sh -q merge_annotation_unique_proteins.fasta -t ${task.cpus} -d \$protein_paths -o \$PWD
         echo "[\$DATE] Executing: /scripts/Aegis3.py --merged_annotation aegis_final_merged_annotations.pkl --hard_masked_genome ${edta_masked_genome} --diamond_hits \${diamond_paths} --intermediate_annotation intermediate_annotation.pkl --final_annotation final_annotation.pkl --export_dir \$PWD --final_export_dir \$PWD --update --source_priority \$protein_names"
         /scripts/Aegis3.py --merged_annotation aegis_final_merged_annotations.pkl --hard_masked_genome ${edta_masked_genome} --diamond_hits \${diamond_paths} --intermediate_annotation intermediate_annotation.pkl --final_annotation final_annotation.pkl --export_dir \$PWD --final_export_dir \$PWD --update --source_priority \$protein_names
        """
}
