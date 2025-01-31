// Define the output directory for intermediate files
// final results are stored in ${projectDir}/FINAL_OUTPUT
params.outdir = "${projectDir}/intermediate_files"

// Define a default log file if not set by the user in nextflow.config
params.logfile = params.logfile ?: "${projectDir}/pipeline_execution.log"

// Print the logfile location at the start of the pipeline
log.info "Logging pipeline execution to: ${params.logfile}"

// Include various processing modules
include { prepare_RNAseq_fastq_files_short } from "./modules/prepare_RNAseq_fastq_files_short"
include { prepare_RNAseq_fastq_files_long } from "./modules/prepare_RNAseq_fastq_files_long"
include { trimming_fastq } from "./modules/trimming_fastq"
include { liftoff_annotations } from "./modules/liftoff_annotations"
include { gffread_convert_gff3_to_cds_fasta } from "./modules/gffread_convert_gff3_to_cds_fasta"
include { salmon_index } from "./modules/salmon_index"
include { salmon_strand_inference } from "./modules/salmon_strand_inference"
include { star_genome_indices } from "./modules/star_genome_indices"
include { star_alignment } from "./modules/star_alignment"
include { hisat2_genome_indices } from "./modules/hisat2_genome_indices"
include { hisat2_alignment } from "./modules/hisat2_alignment"
include { minimap2_genome_indices } from "./modules/minimap2_genome_indices"
include { minimap2_alignment } from "./modules/minimap2_alignment"
include { assembly_transcriptome_star_psiclass } from "./modules/assembly_transcriptome_star_psiclass"
include { assembly_transcriptome_star_stringtie } from "./modules/assembly_transcriptome_star_stringtie"
include { assembly_transcriptome_hisat2_stringtie } from "./modules/assembly_transcriptome_hisat2_stringtie"
include { gffcompare } from "./modules/gffcompare"
include { assembly_transcriptome_minimap2_stringtie } from "./modules/assembly_transcriptome_minimap2_stringtie"
include { Stringtie_merging_short_reads_STAR } from "./modules/Stringtie_merging_short_reads_STAR"
include { Stringtie_merging_short_reads_hisat2 } from "./modules/Stringtie_merging_short_reads_hisat2"
include { Stringtie_merging_long_reads } from "./modules/Stringtie_merging_long_reads"
include { EDTA } from "./modules/EDTA"
include { braker3_prediction } from "./modules/braker3_prediction"
include { braker3_prediction_with_long_reads } from "./modules/braker3_prediction_with_long_reads"
// include { diamond2go } from "./modules/diamond2go"

// Parse RNAseq samplesheet for different types of reads (long, single and paired)
Channel.fromPath( file(params.RNAseq_samplesheet) )
                    .splitCsv(header: true, sep: ',')
                    .filter( ~/.*long.*/ )
                    .map { row ->
                          def sample_ID = row.sampleID
                          def SRA_or_FASTQ = row.SRA_or_FASTQ
                          def library_layout = row.library_layout

                          return [sample_ID, SRA_or_FASTQ, library_layout]
                    }
                    .set{ samples_list_long_reads }

Channel.fromPath( file(params.RNAseq_samplesheet) )
                    .splitCsv(header: true, sep: ',')
                    .filter( ~/.*single.*/ )
                    .map { row ->
                          def sample_ID = row.sampleID
                          def SRA_or_FASTQ = row.SRA_or_FASTQ
                          def library_layout = row.library_layout

                          return [sample_ID, SRA_or_FASTQ, library_layout]
                    }
                    .set{ samples_list_single_short_reads }

Channel.fromPath( file(params.RNAseq_samplesheet) )
                    .splitCsv(header: true, sep: ',')
                    .filter( ~/.*paired.*/ )
                    .map { row ->
                          def sample_ID = row.sampleID
                          def SRA_or_FASTQ = row.SRA_or_FASTQ
                          def library_layout = row.library_layout

                          return [sample_ID, SRA_or_FASTQ, library_layout]
                    }
                    .set{ samples_list_paired_short_reads }

// Combine single and paired short reads into one channel
samples_list_single_short_reads
    .concat(samples_list_paired_short_reads)
    .set { samples_list_short_reads }

// Parse protein samplesheet for BRAKER3
Channel.fromPath( file(params.protein_samplesheet) )
                    .splitCsv(header: true, sep: ',')
                    .map { row ->
                          def organism = row.organism
                          def filename = row.filename

                          return [organism, filename]
                    }
                    .set{ protein_list }

// Workflow definition
workflow {

  // ----------------------------------------------------------------------------------------
  //                           Download/prepare RNAseq reads - OPTIONAL for long reads
  // ----------------------------------------------------------------------------------------

  // Prepare RNAseq short reads for processing
  prepare_RNAseq_fastq_files_short(samples_list_short_reads)

  // Prepare long reads (if any) for processing
  prepare_RNAseq_fastq_files_long(samples_list_long_reads)

  // Trim Illumina short reads
  trimming_fastq(prepare_RNAseq_fastq_files_short.out)

  // ----------------------------------------------------------------------------------------
  //                                Liftoff previous annotations
  // ----------------------------------------------------------------------------------------

  // Lift over previous annotations to new assembly
  liftoff_annotations(file(params.new_assembly).getParent(),file(params.new_assembly).getName(),file(params.previous_assembly).getParent(),file(params.previous_assembly).getName(),file(params.previous_annotations).getParent(),file(params.previous_annotations).getName())

  // -----------------------------------------------------------------------------------------------------------------
  //                                gffread to convert liftoff.gff3 to cds.fasta for Salmon strand inference
  // -----------------------------------------------------------------------------------------------------------------

  // Convert GFF3 to CDS FASTA for Salmon strand inference
  gffread_convert_gff3_to_cds_fasta(file(params.new_assembly).getParent(),file(params.new_assembly).getName(),liftoff_annotations.out.liftoff_previous_annotations)

  // -----------------------------------------------------------------------------------------------------------------------------------------------
  //         Run Salmon for strand inference and classify samples in three strand types : unstranded, stranded_forward and stranded_reverse
  // -----------------------------------------------------------------------------------------------------------------------------------------------
  
  // Salmon index creation and strand inference
  salmon_index(gffread_convert_gff3_to_cds_fasta.out)
  salmon_strand_inference(trimming_fastq.out, salmon_index.out)

  def salmon_output_processed = salmon_strand_inference.out.map { sample_ID, library_layout, reads, strand_file ->
      def strand_info = file(strand_file).text.trim()
      return [sample_ID, library_layout, reads, strand_info]
  }

  // ----------------------------------------------------------------------------------------
  //                    Illumina short RNAseq reads alignment with STAR
  // ----------------------------------------------------------------------------------------

  // Align short RNAseq reads using STAR and strand information + PE/SE information
  star_genome_indices(file(params.new_assembly).getParent(),file(params.new_assembly).getName())

  star_alignment(star_genome_indices.out, salmon_output_processed)

  star_alignment
    .out
    .collect()
    .map { it[0] }
    .set{ concat_star_bams_BRAKER3 }

  // ----------------------------------------------------------------------------------------
  //                    Illumina short RNAseq reads alignment with HISAT2
  // ----------------------------------------------------------------------------------------

  // Align short RNAseq reads using HISAT2 and strand information + PE/SE information
  hisat2_genome_indices(file(params.new_assembly).getParent(),file(params.new_assembly).getName())

  hisat2_alignment(hisat2_genome_indices.out, salmon_output_processed,file(params.new_assembly).getName())

  // ----------------------------------------------------------------------------------------
  //               Pacbio/Nanopore long RNAseq reads alignment with Minimap2  - OPTIONAL
  // ----------------------------------------------------------------------------------------

  // Align long RNAseq reads with Minimap2 (if available)
  minimap2_genome_indices(file(params.new_assembly).getParent(), file(params.new_assembly).getName())
  minimap2_alignment(minimap2_genome_indices.out, prepare_RNAseq_fastq_files_long.out)

  minimap2_alignment
    .out
    .collect()
    .map { it[0] }
    .set { concat_minimap2_bams }

  // ----------------------------------------------------------------------------------------
  //              transcriptome assembly with PsiCLASS on STAR alignments (short reads)
  // ----------------------------------------------------------------------------------------

  // Transcriptome assembly with PsiCLASS
  assembly_transcriptome_star_psiclass(star_alignment.out) 

  // ----------------------------------------------------------------------------------------
  //        transcriptome assembly with Stringtie on STAR alignments (short reads)
  // ----------------------------------------------------------------------------------------
  // Transcriptome assembly with Stringtie

  assembly_transcriptome_star_stringtie(star_alignment.out)

  assembly_transcriptome_star_stringtie.out.star_stringtie_transcriptome_gtf
  .groupTuple()
  .collect()
  .map { it[0] }
  .set { concat_star_stringtie_for_merging }

  // Stringtie merging of all short reads transcriptomes (STAR/Stringtie)
  Stringtie_merging_short_reads_STAR(concat_star_stringtie_for_merging)

  // ----------------------------------------------------------------------------------------
  //        transcriptome assembly with Stringtie on HISAT2 alignments (short reads)
  // ----------------------------------------------------------------------------------------
  // Transcriptome assembly with Stringtie
  
  assembly_transcriptome_hisat2_stringtie(hisat2_alignment.out)

  assembly_transcriptome_hisat2_stringtie.out.hisat2_stringtie_transcriptome_gtf
  .groupTuple()
  .collect()
  .map { it[0] }
  .set { concat_hisat2_stringtie_for_merging }

  // Stringtie merging of all short reads transcriptomes (HISAT2/Stringtie)
  Stringtie_merging_short_reads_hisat2(concat_hisat2_stringtie_for_merging)

  // ----------------------------------------------------------------------------------------
  //        gffcompare to merge PsiCLASS transcriptomes
  // ----------------------------------------------------------------------------------------

  assembly_transcriptome_star_psiclass.out
    .collect()
    .map { it[0] }
    .set { concat_star_psiclass_for_merging }

  // GFFcompare to merge PsiCLASS transcriptomes
  gffcompare(concat_star_psiclass_for_merging)

  // ----------------------------------------------------------------------------------------
  //      transcriptome assembly with Stringtie on minimap2 alignments (long reads) - OPTIONAL
  // ----------------------------------------------------------------------------------------

  // Transcriptome assembly with StringTie on Minimap2 alignments (if long reads)
  assembly_transcriptome_minimap2_stringtie(minimap2_alignment.out)

  assembly_transcriptome_minimap2_stringtie.out.minimap2_stringtie_transcriptome_gtf
  .groupTuple()
  .collect()
  .map { it[0] }
  .set { concat_minimap2_stringtie_for_merging }

  Stringtie_merging_long_reads(concat_minimap2_stringtie_for_merging)

  // ----------------------------------------------------------------------------------------
  // -------------------------- Genome masking with EDTA ------------------------------------
  // ----------------------------------------------------------------------------------------

  // Optionally run EDTA for genome masking
  if (params.EDTA == 'yes') {
    println "Running EDTA process"
    EDTA(file(params.new_assembly).getParent(), file(params.new_assembly).getName())
  } else {
    println "Skipping EDTA process (params.EDTA = '${params.EDTA}'). To launch EDTA, put EDTA = 'yes' in nextflow.config file."
  }

  // ----------------------------------------------------------------------------------------
  //                                    BRAKER3 (AUGUSTUS/Genemark)
  // ----------------------------------------------------------------------------------------

  // Gene prediction using BRAKER3 with or without long reads

  braker3_prediction(
      file(params.new_assembly).getParent(),
      file(params.new_assembly).getName(),
      file(params.protein_samplesheet).getParent(),
      file(params.protein_samplesheet).getName(),
      concat_star_bams_BRAKER3
  )

  braker3_prediction_with_long_reads(
    file(params.new_assembly).getParent(),
    file(params.new_assembly).getName(),
    file(params.protein_samplesheet).getParent(),
    file(params.protein_samplesheet).getName(),
    concat_star_bams_BRAKER3,
    concat_minimap2_bams
  )

  // ----------------------------------------------------------------------------------------
  //     Aegis scripts (1, 2, 3) to create the final GFF3 file from all the evidences
  // ----------------------------------------------------------------------------------------

    // Placeholder for Aegis scripts to create the final GFF3 file
    // TO DO

  // ----------------------------------------------------------------------------------------
  //               Diamond2GO on proteins predicted with TITAN
  // ----------------------------------------------------------------------------------------
  // diamond2go(proteins_file)
}
