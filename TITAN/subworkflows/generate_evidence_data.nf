nextflow.enable.dsl = 2

// Include various processing modules
include { prepare_RNAseq_fastq_files_short } from "../modules/prepare_RNAseq_fastq_files_short"
include { prepare_RNAseq_fastq_files_long } from "../modules/prepare_RNAseq_fastq_files_long"
include { trimming_fastq } from "../modules/trimming_fastq"
include { liftoff_annotations } from "../modules/liftoff_annotations"
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

    main:
        
        // Prepare RNAseq short reads for processing
        short_reads_prepared = prepare_RNAseq_fastq_files_short(samples_list_short_reads)
        
        // Prepare long reads (if any) for processing
        long_reads_prepared = prepare_RNAseq_fastq_files_long(samples_list_long_reads)
        
        // Trim Illumina short reads
        trimmed_reads = trimming_fastq(short_reads_prepared.prepared_fastqs)
        
        // Lift over previous annotations to new assembly
        previous_annotations = liftoff_annotations(
            file(params.new_assembly).getParent(),
            file(params.new_assembly).getName(),
            file(params.previous_assembly).getParent(),
            file(params.previous_assembly).getName(),
            file(params.previous_annotations).getParent(),
            file(params.previous_annotations).getName()
        )
        
        // Convert GFF3 to CDS FASTA for Salmon strand inference
        gff_cds = agat_convert_gff3_to_cds_fasta(
            file(params.new_assembly).getParent(),
            file(params.new_assembly).getName(),
            previous_annotations.liftoff_previous_annotations
        )
        
        // Salmon index and strand inference
        salmon_index_result = salmon_index(gff_cds)
        strand_inference = salmon_strand_inference(trimmed_reads.trimmed_reads, salmon_index_result.index)
        
        // Process strand inference results
        salmon_output_processed = strand_inference.strand_inference_result.map { sample_ID, library_layout, reads, strand_file ->
            def strand_info = file(strand_file).text.trim()
            return [sample_ID, library_layout, reads, strand_info]
        }
        
        // Process STAR alignments
        star_indices = star_genome_indices(
            file(params.new_assembly).getParent(),
            file(params.new_assembly).getName()
        )
        star_aligned = star_alignment(trimmed_reads.trimmed_reads, star_indices.index)
        
        // Process HISAT2 alignments
        hisat2_indices = hisat2_genome_indices(
            file(params.new_assembly).getParent(),
            file(params.new_assembly).getName()
        )
        hisat2_aligned = hisat2_alignment(hisat2_indices.index, salmon_output_processed, file(params.new_assembly).getName())
        
        // Process hisat2 assemblies
        hisat2_assemblies_stringtie = assembly_transcriptome_hisat2_stringtie(hisat2_aligned.samples_aligned)

        concat_hisat2_stringtie_for_merging = hisat2_assemblies_stringtie.hisat2_stringtie_transcriptome_gtf
        .groupTuple()
        .collect()
        .map { it[0] }

        merged_hisat2_stringtie = Stringtie_merging_short_reads_hisat2(concat_hisat2_stringtie_for_merging)

        // Process minimap2 alignments for long reads
        if (params.use_long_reads) {
            minimap2_indices = minimap2_genome_indices(
                file(params.new_assembly).getParent(),
                file(params.new_assembly).getName()
            )
            minimap2_aligned = minimap2_alignment(long_reads_prepared.prepared_fastqs, minimap2_indices.index)
            
            // Process long read assemblies
            long_reads_assemblies = assembly_transcriptome_minimap2_stringtie(minimap2_aligned.samples_aligned)

            concat_minimap2_stringtie_for_merging = long_reads_assemblies.minimap2_stringtie_transcriptome_gtf
            .groupTuple()
            .collect()
            .map { it[0] }

            merged_long_reads = Stringtie_merging_long_reads(concat_minimap2_stringtie_for_merging)
        }
        
        // Process short read assemblies
        star_assemblies_stringtie = assembly_transcriptome_star_stringtie(star_aligned.samples_aligned)
        star_assemblies_psiclass = assembly_transcriptome_star_psiclass(star_aligned.samples_aligned)
        
        concat_star_stringtie_for_merging = star_assemblies_stringtie.star_stringtie_transcriptome_gtf
        .groupTuple()
        .collect()
        .map { it[0] }

        // Merge assemblies
        merged_star_stringtie = Stringtie_merging_short_reads_STAR(concat_star_stringtie_for_merging)

        concat_star_psiclass_for_merging = star_assemblies_psiclass.psiclass_assemblies
          .collect()
          .map { it[0] }

        // GFFcompare to merge PsiCLASS transcriptomes
        gffcompare_out = gffcompare(concat_star_psiclass_for_merging)

        // Run EDTA if enabled
        if (params.EDTA == 'yes') {
            masked_genome = EDTA(
                file(params.new_assembly).getParent(),
                file(params.new_assembly).getName()
            )
        }
        
        // Run BRAKER3
        if (params.use_long_reads) {
            braker3_results = braker3_prediction_with_long_reads(
                file(params.new_assembly).getParent(),
                file(params.new_assembly).getName(),
                file(params.protein_samplesheet).getParent(),
                file(params.protein_samplesheet).getName(),
                star_aligned.samples_aligned,
                minimap2_aligned.samples_aligned
            )
        } else {
            braker3_results = braker3_prediction(
                file(params.new_assembly).getParent(),
                file(params.new_assembly).getName(),
                file(params.protein_samplesheet).getParent(),
                file(params.protein_samplesheet).getName(),
                star_aligned.samples_aligned
            )
        }

    // emit outputs
    outputs_ch = Channel.empty()

    // Add EDTA outputs if enabled
    if (params.EDTA == 'yes') {
        outputs_ch = outputs_ch.mix(
            Channel.of("masked_genome.masked_genome").combine(masked_genome.masked_genome)
        )
    }

    // Add long reads outputs if enabled
    if (params.use_long_reads) {
        outputs_ch = outputs_ch.mix(
            Channel.of("merged_long_reads.default_args_gff").combine(merged_long_reads.default_args_gff),
            Channel.of("merged_long_reads.alt_args_gff").combine(merged_long_reads.alt_args_gff)
        )
    }

    outputs_ch = outputs_ch.mix(
        Channel.of("braker3_results.augustus_gff").combine(braker3_results.augustus_gff),
        Channel.of("braker3_results.genemark_gtf").combine(braker3_results.genemark_gtf),
        Channel.of("previous_annotations.liftoff_previous_annotations").combine(previous_annotations.liftoff_previous_annotations),
        Channel.of("merged_star_stringtie.default_args_stranded").combine(merged_star_stringtie.default_args_stranded),
        Channel.of("merged_star_stringtie.alt_args_stranded").combine(merged_star_stringtie.alt_args_stranded),
        Channel.of("gffcompare_out.star_psiclass_stranded").combine(gffcompare_out.star_psiclass_stranded),
        Channel.of("merged_star_stringtie.default_args_unstranded").combine(merged_star_stringtie.default_args_unstranded),
        Channel.of("merged_star_stringtie.alt_args_unstranded").combine(merged_star_stringtie.alt_args_unstranded),
        Channel.of("gffcompare_out.star_psiclass_unstranded").combine(gffcompare_out.star_psiclass_unstranded)
    )

    emit:
        results = outputs_ch
}