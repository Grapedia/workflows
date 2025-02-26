nextflow.enable.dsl=2

def output_directory = params.output_dir

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

    // Define a dictionnary (Map) to store outputs for aegis workflow
    def outputs_map = [:]
    
    // ----------------------------------------------------------------------------------------
    //                           Download/prepare RNAseq reads - OPTIONAL for long reads
    // ----------------------------------------------------------------------------------------

    // Prepare RNAseq short reads for processing
    short_reads_prepared = prepare_RNAseq_fastq_files_short(samples_list_short_reads)

    // Prepare long reads (if any) for processing
    long_reads_prepared = prepare_RNAseq_fastq_files_long(samples_list_long_reads)

    // Trim Illumina short reads
    trimmed_reads = trimming_fastq(short_reads_prepared.prepared_fastqs)

    // ----------------------------------------------------------------------------------------
    //                                Liftoff previous annotations
    // ----------------------------------------------------------------------------------------

    // Lift over previous annotations to new assembly
    previous_annotations = liftoff_annotations(file(params.new_assembly).getParent(),file(params.new_assembly).getName(),file(params.previous_assembly).getParent(),file(params.previous_assembly).getName(),file(params.previous_annotations).getParent(),file(params.previous_annotations).getName()).collect()
    println "Value of previous_annotations : ${previous_annotations}"

    // -----------------------------------------------------------------------------------------------------------------
    //                                gffread to convert liftoff.gff3 to cds.fasta for Salmon strand inference
    // -----------------------------------------------------------------------------------------------------------------

    // Convert GFF3 to CDS FASTA for Salmon strand inference
    gff_cds = agat_convert_gff3_to_cds_fasta(file(params.new_assembly).getParent(),file(params.new_assembly).getName(),previous_annotations.liftoff_previous_annotations)

    // -----------------------------------------------------------------------------------------------------------------------------------------------
    //         Run Salmon for strand inference and classify samples in three strand types : unstranded, stranded_forward and stranded_reverse
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    
    // Salmon index creation and strand inference
    salmon_index_result = salmon_index(gff_cds.gffread_cds)
    strand_inference = salmon_strand_inference(trimmed_reads.trimmed_reads, salmon_index_result.index)

    salmon_output_processed = strand_inference.strand_inference_result.map { sample_ID, library_layout, reads, strand_file ->
        def strand_info = file(strand_file).text.trim()
        return [sample_ID, library_layout, reads, strand_info]
    }

    // ----------------------------------------------------------------------------------------
    //                    Illumina short RNAseq reads alignment with STAR
    // ----------------------------------------------------------------------------------------

    // Align short RNAseq reads using STAR and strand information + PE/SE information
    star_indices = star_genome_indices(file(params.new_assembly).getParent(),file(params.new_assembly).getName())

    star_aligned = star_alignment(star_indices.out, salmon_output_processed)

    concat_star_bams_BRAKER3 = star_aligned
      .samples_aligned
      .collect()
      .map { it[0] }

    // ----------------------------------------------------------------------------------------
    //                    Illumina short RNAseq reads alignment with HISAT2
    // ----------------------------------------------------------------------------------------

    // Align short RNAseq reads using HISAT2 and strand information + PE/SE information
    hisat2_indices = hisat2_genome_indices(file(params.new_assembly).getParent(),file(params.new_assembly).getName())

    hisat2_aligned = hisat2_alignment(hisat2_indices.index, salmon_output_processed,file(params.new_assembly).getName())

    // ----------------------------------------------------------------------------------------
    //               Pacbio/Nanopore long RNAseq reads alignment with Minimap2  - OPTIONAL
    // ----------------------------------------------------------------------------------------

    // Align long RNAseq reads with Minimap2 (if available)
    minimap2_indices = minimap2_genome_indices(file(params.new_assembly).getParent(), file(params.new_assembly).getName())
    minimap2_aligned = minimap2_alignment(minimap2_indices.index, long_reads_prepared.prepared_samples)

    concat_minimap2_bams = minimap2_aligned
      .samples_aligned
      .collect()
      .map { it[0] }

    // ----------------------------------------------------------------------------------------
    //              transcriptome assembly with PsiCLASS on STAR alignments (short reads)
    // ----------------------------------------------------------------------------------------

    // Transcriptome assembly with PsiCLASS
    transcriptome_star_psiclass = assembly_transcriptome_star_psiclass(star_aligned.samples_aligned) 

    // ----------------------------------------------------------------------------------------
    //        transcriptome assembly with Stringtie on STAR alignments (short reads)
    // ----------------------------------------------------------------------------------------
    // Transcriptome assembly with Stringtie

    transcriptome_star_stringtie = assembly_transcriptome_star_stringtie(star_aligned.samples_aligned)

    concat_star_stringtie_for_merging = transcriptome_star_stringtie.star_stringtie_transcriptome_gtf
    .groupTuple()
    .collect()
    .map { it[0] }

    // Stringtie merging of all short reads transcriptomes (STAR/Stringtie)
    Stringtie_merging_short_reads_STAR_out = Stringtie_merging_short_reads_STAR(concat_star_stringtie_for_merging).collect()
    println "Value of Stringtie_merging_short_reads_STAR_out : ${Stringtie_merging_short_reads_STAR_out}"

    // ----------------------------------------------------------------------------------------
    //        transcriptome assembly with Stringtie on HISAT2 alignments (short reads)
    // ----------------------------------------------------------------------------------------
    // Transcriptome assembly with Stringtie
    
    transcriptome_hisat2_stringtie = assembly_transcriptome_hisat2_stringtie(hisat2_aligned.samples_aligned)

    concat_hisat2_stringtie_for_merging = transcriptome_hisat2_stringtie.hisat2_stringtie_transcriptome_gtf
    .groupTuple()
    .collect()
    .map { it[0] }

    // Stringtie merging of all short reads transcriptomes (HISAT2/Stringtie)
    Stringtie_merging_short_reads_hisat2_out = Stringtie_merging_short_reads_hisat2(concat_hisat2_stringtie_for_merging)

    // ----------------------------------------------------------------------------------------
    //        gffcompare to merge PsiCLASS transcriptomes
    // ----------------------------------------------------------------------------------------

    concat_star_psiclass_for_merging = transcriptome_star_psiclass.psiclass_assemblies
      .collect()
      .map { it[0] }

    // GFFcompare to merge PsiCLASS transcriptomes
    gffcompare_out = gffcompare(concat_star_psiclass_for_merging).collect()
    println "Value of gffcompare_out : ${gffcompare_out}"

    // ----------------------------------------------------------------------------------------
    //      transcriptome assembly with Stringtie on minimap2 alignments (long reads) - OPTIONAL
    // ----------------------------------------------------------------------------------------

    // Transcriptome assembly with StringTie on Minimap2 alignments (if long reads)
    transcriptome_minimap2_stringtie = assembly_transcriptome_minimap2_stringtie(minimap2_aligned.samples_aligned)

    concat_minimap2_stringtie_for_merging = transcriptome_minimap2_stringtie.minimap2_stringtie_transcriptome_gtf
    .groupTuple()
    .collect()
    .map { it[0] }

    Stringtie_merging_long_reads_out = Stringtie_merging_long_reads(concat_minimap2_stringtie_for_merging).collect()
    println "Value of Stringtie_merging_long_reads_out : ${Stringtie_merging_long_reads_out}"

    // ----------------------------------------------------------------------------------------
    // -------------------------- Genome masking with EDTA ------------------------------------
    // ----------------------------------------------------------------------------------------

    // Optionally run EDTA for genome masking
    if (params.EDTA == 'yes') {
      println "Running EDTA process"
      masked_genome = EDTA(file(params.new_assembly).getParent(), file(params.new_assembly).getName()).collect()
      println "Value of masked_genome : ${masked_genome}"
    } else {
      println "Skipping the EDTA and AEGIS processes because EDTA was not run (params.EDTA = '${params.EDTA}'). To enable EDTA, and consequently AEGIS, set EDTA = 'yes' in the nextflow.config file."
    }

    // ----------------------------------------------------------------------------------------
    //                                    BRAKER3 (AUGUSTUS/Genemark)
    // ----------------------------------------------------------------------------------------

    // Gene prediction using BRAKER3 with or without long reads

    if (params.use_long_reads) {

      braker3_results = braker3_prediction_with_long_reads(
        file(params.new_assembly).getParent(),
        file(params.new_assembly).getName(),
        file(params.protein_samplesheet).getParent(),
        file(params.protein_samplesheet).getName(),
        concat_star_bams_BRAKER3,
        concat_minimap2_bams
      ).collect()
      println "Value of braker3_results : ${braker3_results}"

    } else {

      braker3_results = braker3_prediction(
          file(params.new_assembly).getParent(),
          file(params.new_assembly).getName(),
          file(params.protein_samplesheet).getParent(),
          file(params.protein_samplesheet).getName(),
          concat_star_bams_BRAKER3
      ).collect()
      println "Value of braker3_results : ${braker3_results}"

    }

    // masked genome
    if (params.EDTA == 'yes') {
      outputs_map["masked_genome"] = masked_genome.masked_genome
    }

    // BRAKER3 results
    outputs_map["braker3_augustus"] = braker3_results.augustus_gff
    outputs_map["braker3_genemark"] = braker3_results.genemark_gff

    // Outputs from long reads if "true"
    if (params.use_long_reads) {
      outputs_map["stringtie_long_reads_default"] = Stringtie_merging_long_reads_out.default_args_gff
      outputs_map["stringtie_long_reads_alt"] = Stringtie_merging_long_reads_out.alt_args_gff
    }

    // Liftoff annotations
    outputs_map["liftoff_annotations"] = previous_annotations.liftoff_previous_annotations

    // Results of StringTie on short reads STAR for stranded and unstranded data
    outputs_map["stringtie_short_reads_STAR_default_stranded"] = Stringtie_merging_short_reads_STAR_out.default_args_stranded
    outputs_map["stringtie_short_reads_STAR_alt_stranded"] = Stringtie_merging_short_reads_STAR_out.alt_args_stranded
    outputs_map["stringtie_short_reads_STAR_default_unstranded"] = Stringtie_merging_short_reads_STAR_out.default_args_unstranded
    outputs_map["stringtie_short_reads_STAR_alt_unstranded"] = Stringtie_merging_short_reads_STAR_out.alt_args_unstranded

    // Results of gffcompare after PsiClass for stranded and unstranded data
    outputs_map["gffcompare_star_psiclass_stranded"] = gffcompare_out.star_psiclass_stranded
    outputs_map["gffcompare_star_psiclass_unstranded"] = gffcompare_out.star_psiclass_unstranded

    println "Debug - outputs_map content:"
    println outputs_map

  emit:
    outputs_map
}