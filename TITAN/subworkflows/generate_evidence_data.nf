nextflow.enable.dsl = 2

// Include various processing modules
include { prepare_RNAseq_fastq_files_short } from "../modules/prepare_RNAseq_fastq_files_short"
include { prepare_RNAseq_fastq_files_long } from "../modules/prepare_RNAseq_fastq_files_long"
include { trimming_fastq } from "../modules/trimming_fastq"
include { liftoff_annotations } from "../modules/liftoff_annotations"
include { egapx } from "../modules/egapx"
include { agat_convert_gff3_to_cds_fasta } from "../modules/agat_convert_gff3_to_cds_fasta"
include { salmon_index } from "../modules/salmon_index"
include { salmon_strand_inference } from "../modules/salmon_strand_inference"
include { star_genome_indices } from "../modules/star_genome_indices"
include { star_alignment } from "../modules/star_alignment"
include { hisat2_genome_indices } from "../modules/hisat2_genome_indices"
include { hisat2_alignment } from "../modules/hisat2_alignment"
include { minimap2_genome_indices } from "../modules/minimap2_genome_indices"
include { minimap2_alignment } from "../modules/minimap2_alignment"
include { assembly_transcriptome_star_psiclass } from "../modules/assembly_transcriptome_star_psiclass"
include { assembly_transcriptome_star_stringtie } from "../modules/assembly_transcriptome_star_stringtie"
include { assembly_transcriptome_hisat2_stringtie } from "../modules/assembly_transcriptome_hisat2_stringtie"
include { gffcompare } from "../modules/gffcompare"
include { assembly_transcriptome_minimap2_stringtie } from "../modules/assembly_transcriptome_minimap2_stringtie"
include { Stringtie_merging_short_reads_STAR } from "../modules/Stringtie_merging_short_reads_STAR"
include { Stringtie_merging_short_reads_hisat2 } from "../modules/Stringtie_merging_short_reads_hisat2"
include { Stringtie_merging_long_reads } from "../modules/Stringtie_merging_long_reads"
include { EDTA } from "../modules/EDTA"
include { braker3_prediction } from "../modules/braker3_prediction"
include { braker3_prediction_with_long_reads } from "../modules/braker3_prediction_with_long_reads"

workflow generate_evidence_data {
    take:
        samples_list_long_reads
        samples_list_short_reads
        protein_list
        has_long_reads

    main:
        def new_assembly = file(params.new_assembly)
        def previous_assembly = file(params.previous_assembly)
        def previous_annotations_file = file(params.previous_annotations)
        def edta_script = file("${projectDir}/scripts/edta.sh")
        def stringtie_script = file("${projectDir}/scripts/Stringtie.sh")
        def stringtie_alt_script = file("${projectDir}/scripts/Stringtie_AltCommands.sh")
        def stringtie_transcriptome_script = file("${projectDir}/scripts/run_stringtie_transcriptome.sh")
        def empty_gtf = file("${projectDir}/assets/empty.gtf")
        def empty_default_gtf = file("${projectDir}/assets/empty_default.gtf")
        def empty_alt_gtf = file("${projectDir}/assets/empty_alt.gtf")
        def empty_psiclass_gtf = file("${projectDir}/assets/empty_psiclass.gtf")

        // Prepare RNAseq short reads for processing
        short_reads_prepared = prepare_RNAseq_fastq_files_short(
            samples_list_short_reads,
            params.ena_download_timeout_seconds,
            params.ena_max_download_attempts,
            params.ena_retry_wait_seconds,
            params.ena_verify_md5
        )

        // Trim Illumina short reads
        trimmed_reads = trimming_fastq(short_reads_prepared.prepared_fastqs)

        // Lift over previous annotations to new assembly
        previous_annotations = liftoff_annotations(
            new_assembly,
            new_assembly.getName(),
            previous_assembly,
            previous_assembly.getName(),
            previous_annotations_file,
            previous_annotations_file.getName()
        )

        // EGAPx annotation pipeline on new assembly
        egapx_annotations = egapx(file(params.egapx_paramfile))

        // Convert GFF3 to CDS FASTA for Salmon strand inference
        gff_cds = agat_convert_gff3_to_cds_fasta(
            file(params.new_assembly),
            previous_annotations.liftoff_previous_annotations
        )

        // Salmon index and strand inference
        salmon_index_result = salmon_index(gff_cds.cds_fasta)
        strand_inference = salmon_strand_inference(trimmed_reads.trimmed_reads, salmon_index_result.index)

        // Process strand inference results
        salmon_output_processed = strand_inference.strand_inference_result.map { sample_ID, library_layout, reads, strand_info, strand_file ->
            return [sample_ID, library_layout, reads, strand_info]
        }

        // Process STAR alignments
        star_indices = star_genome_indices(
            new_assembly,
            new_assembly.getName()
        )
        star_aligned = star_alignment(star_indices.index, salmon_output_processed, params.STAR_memory_per_job)

        star_aligned.samples_aligned
          .map { sample_ID, bam_file, strand_type -> bam_file }
          .collect()
          .set{ concat_star_bams_BRAKER3 }

        // Process HISAT2 alignments
        hisat2_indices = hisat2_genome_indices(
            new_assembly,
            new_assembly.getName()
        )
        hisat2_aligned = hisat2_alignment(hisat2_indices.index, salmon_output_processed, file(params.new_assembly).getName())

        // Process hisat2 assemblies
        hisat2_assemblies_stringtie = assembly_transcriptome_hisat2_stringtie(
            hisat2_aligned.samples_aligned,
            stringtie_script,
            stringtie_alt_script,
            stringtie_transcriptome_script
        )

        hisat2_assemblies_stringtie.hisat2_stringtie_transcriptomes
          .filter { sample_ID, default_gtf, alt_gtf, strand_type -> strand_type != 'unstranded' }
          .map { sample_ID, default_gtf, alt_gtf, strand_type -> default_gtf }
          .collect()
          .set { hisat2_stranded_default_gtfs }

        hisat2_assemblies_stringtie.hisat2_stringtie_transcriptomes
          .filter { sample_ID, default_gtf, alt_gtf, strand_type -> strand_type != 'unstranded' }
          .map { sample_ID, default_gtf, alt_gtf, strand_type -> alt_gtf }
          .collect()
          .set { hisat2_stranded_alt_gtfs }

        hisat2_assemblies_stringtie.hisat2_stringtie_transcriptomes
          .filter { sample_ID, default_gtf, alt_gtf, strand_type -> strand_type == 'unstranded' }
          .map { sample_ID, default_gtf, alt_gtf, strand_type -> default_gtf }
          .collect()
          .ifEmpty(empty_default_gtf)
          .set { hisat2_unstranded_default_gtfs }

        hisat2_assemblies_stringtie.hisat2_stringtie_transcriptomes
          .filter { sample_ID, default_gtf, alt_gtf, strand_type -> strand_type == 'unstranded' }
          .map { sample_ID, default_gtf, alt_gtf, strand_type -> alt_gtf }
          .collect()
          .ifEmpty(empty_alt_gtf)
          .set { hisat2_unstranded_alt_gtfs }

        merged_hisat2_stringtie = Stringtie_merging_short_reads_hisat2(
            hisat2_stranded_default_gtfs,
            hisat2_stranded_alt_gtfs,
            hisat2_unstranded_default_gtfs,
            hisat2_unstranded_alt_gtfs
        )

        merged_long_reads_default_args_gff = Channel.value(empty_default_gtf)
        merged_long_reads_alt_args_gff = Channel.value(empty_alt_gtf)

        // Process minimap2 alignments for long reads
        if (has_long_reads) {
            // Prepare long reads (if any) for processing
            long_reads_prepared = prepare_RNAseq_fastq_files_long(
                samples_list_long_reads,
                params.ena_download_timeout_seconds,
                params.ena_max_download_attempts,
                params.ena_retry_wait_seconds,
                params.ena_verify_md5
            )

            minimap2_indices = minimap2_genome_indices(
                new_assembly,
                new_assembly.getName()
            )
            minimap2_aligned = minimap2_alignment(minimap2_indices.index, long_reads_prepared.prepared_fastqs)

            minimap2_aligned.samples_aligned
              .map { sample_ID, bam_file -> bam_file }
              .collect()
              .set { concat_minimap2_bams_BRAKER3 }

            // Process long read assemblies
            long_reads_assemblies = assembly_transcriptome_minimap2_stringtie(
                minimap2_aligned.samples_aligned,
                stringtie_script,
                stringtie_alt_script,
                stringtie_transcriptome_script
            )

            long_reads_assemblies.minimap2_stringtie_transcriptomes
              .map { sample_ID, default_gtf, alt_gtf -> default_gtf }
              .collect()
              .set { minimap2_default_gtfs }

            long_reads_assemblies.minimap2_stringtie_transcriptomes
              .map { sample_ID, default_gtf, alt_gtf -> alt_gtf }
              .collect()
              .set { minimap2_alt_gtfs }

            merged_long_reads = Stringtie_merging_long_reads(minimap2_default_gtfs, minimap2_alt_gtfs)
            merged_long_reads_default_args_gff = merged_long_reads.default_args_gff
            merged_long_reads_alt_args_gff = merged_long_reads.alt_args_gff
        }

        // Process short read assemblies
        star_assemblies_stringtie = assembly_transcriptome_star_stringtie(
            star_aligned.samples_aligned,
            stringtie_script,
            stringtie_alt_script,
            stringtie_transcriptome_script
        )
        star_assemblies_psiclass = assembly_transcriptome_star_psiclass(
            star_aligned.samples_aligned,
            params.PSICLASS_vd_option,
            params.PSICLASS_c_option
        )

        star_assemblies_stringtie.star_stringtie_transcriptomes
          .filter { sample_ID, default_gtf, alt_gtf, strand_type -> strand_type != 'unstranded' }
          .map { sample_ID, default_gtf, alt_gtf, strand_type -> default_gtf }
          .collect()
          .set { star_stranded_default_gtfs }

        star_assemblies_stringtie.star_stringtie_transcriptomes
          .filter { sample_ID, default_gtf, alt_gtf, strand_type -> strand_type != 'unstranded' }
          .map { sample_ID, default_gtf, alt_gtf, strand_type -> alt_gtf }
          .collect()
          .set { star_stranded_alt_gtfs }

        star_assemblies_stringtie.star_stringtie_transcriptomes
          .filter { sample_ID, default_gtf, alt_gtf, strand_type -> strand_type == 'unstranded' }
          .map { sample_ID, default_gtf, alt_gtf, strand_type -> default_gtf }
          .collect()
          .ifEmpty(empty_default_gtf)
          .set { star_unstranded_default_gtfs }

        star_assemblies_stringtie.star_stringtie_transcriptomes
          .filter { sample_ID, default_gtf, alt_gtf, strand_type -> strand_type == 'unstranded' }
          .map { sample_ID, default_gtf, alt_gtf, strand_type -> alt_gtf }
          .collect()
          .ifEmpty(empty_alt_gtf)
          .set { star_unstranded_alt_gtfs }

        // Merge assemblies
        merged_star_stringtie = Stringtie_merging_short_reads_STAR(
            star_stranded_default_gtfs,
            star_stranded_alt_gtfs,
            star_unstranded_default_gtfs,
            star_unstranded_alt_gtfs
        )

        star_assemblies_psiclass.psiclass_assemblies
          .filter { sample_ID, gtf_file, strand_type -> strand_type != 'unstranded' }
          .map { sample_ID, gtf_file, strand_type -> gtf_file }
          .collect()
          .set { star_psiclass_stranded_gtfs }

        star_assemblies_psiclass.psiclass_assemblies
          .filter { sample_ID, gtf_file, strand_type -> strand_type == 'unstranded' }
          .map { sample_ID, gtf_file, strand_type -> gtf_file }
          .collect()
          .ifEmpty(empty_psiclass_gtf)
          .set { star_psiclass_unstranded_gtfs }

        // GFFcompare to merge PsiCLASS transcriptomes
        gffcompare_out = gffcompare(star_psiclass_stranded_gtfs, star_psiclass_unstranded_gtfs)

        // EDTA is mandatory: Aegis requires the hard-masked genome.
        edta_results = EDTA(
            new_assembly,
            new_assembly.getName(),
            edta_script
        )

        protein_fastas = protein_list
          .map { organism, filename -> file(filename) }
          .collect()

        clean_protein_script = file("${projectDir}/scripts/clean_protein_fasta_for_BRAKER3.py")

        // Run BRAKER3
        if (has_long_reads) {
            braker3_results = braker3_prediction_with_long_reads(
                file(params.new_assembly),
                protein_fastas,
                concat_star_bams_BRAKER3,
                concat_minimap2_bams_BRAKER3,
                clean_protein_script
            )
        } else {
            braker3_results = braker3_prediction(
                file(params.new_assembly),
                protein_fastas,
                concat_star_bams_BRAKER3,
                clean_protein_script
            )
        }

    emit:
        masked_genome = edta_results.masked_genome
        egapx_gff3 = egapx_annotations.gff3
        egapx_gtf = egapx_annotations.gtf
        egapx_proteins = egapx_annotations.proteins
        egapx_cds = egapx_annotations.cds
        egapx_transcripts = egapx_annotations.transcripts
        egapx_annotated_genome_asn = egapx_annotations.annotated_genome_asn
        egapx_output_dir = egapx_annotations.output_dir
        egapx_versions = egapx_annotations.versions
        liftoff_gff3 = previous_annotations.liftoff_previous_annotations
        liftoff_unmapped_features = previous_annotations.unmapped_features
        braker_augustus_gff3 = braker3_results.augustus_gff
        braker_genemark_gtf = braker3_results.genemark_gtf
        braker_genemark_supported_gtf = braker3_results.genemark_supported_gtf
        braker_gff3 = braker3_results.braker_gff
        star_stringtie_stranded_default_gtf = merged_star_stringtie.default_args_stranded
        star_stringtie_stranded_alt_gtf = merged_star_stringtie.alt_args_stranded
        star_stringtie_unstranded_default_gtf = merged_star_stringtie.default_args_unstranded
        star_stringtie_unstranded_alt_gtf = merged_star_stringtie.alt_args_unstranded
        star_psiclass_stranded_gtf = gffcompare_out.star_psiclass_stranded
        star_psiclass_unstranded_gtf = gffcompare_out.star_psiclass_unstranded
        long_reads_default_gtf = merged_long_reads_default_args_gff
        long_reads_alt_gtf = merged_long_reads_alt_args_gff
}
