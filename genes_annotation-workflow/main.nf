// Path to outdir
params.outdir = "${projectDir}/intermediate_files"

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

samples_list_long_reads
    .count()
    .map { count -> count > 0 }
    .set { has_long_reads }

has_long_reads.view { flag ->
    println "Debug: has_long_reads is ${flag}"
}

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

samples_list_single_short_reads.concat(samples_list_paired_short_reads).set { samples_list_short_reads }

Channel.fromPath( file(params.protein_samplesheet) )
                    .splitCsv(header: true, sep: ',')
                    .map { row ->
                          def organism = row.organism
                          def filename = row.filename

                          return [organism, filename]
                    }
                    .set{ protein_list }

workflow {

  // ----------------------------------------------------------------------------------------
  //                           Download/prepare RNAseq reads - OPTIONAL for long reads
  // ----------------------------------------------------------------------------------------
  prepare_RNAseq_fastq_files_short(samples_list_short_reads) // VALIDATED

  // Check that samples_list_long_reads is empty or not before running prepare_RNAseq_fastq_files_long
  prepare_RNAseq_fastq_files_long(samples_list_long_reads)

  has_long_reads.view { flag ->
      if (!flag) {
          println "No long reads samples detected, skipping long reads processing."
      } else {
          println "Long reads samples detected, processing them with prepare_RNAseq_fastq_files_long."
      }
  }

  // trimming with fastp in only done on Illumina short reads
  trimming_fastq(prepare_RNAseq_fastq_files_short.out) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                                Liftoff previous annotations
  // ----------------------------------------------------------------------------------------
  liftoff_annotations(file(params.new_assembly).getParent(),file(params.new_assembly).getName(),file(params.previous_assembly).getParent(),file(params.previous_assembly).getName(),file(params.previous_annotations).getParent(),file(params.previous_annotations).getName()) // VALIDATED

  // -----------------------------------------------------------------------------------------------------------------
  //                                gffread to convert liftoff.gff3 to cds.fasta for Salmon strand inference
  // -----------------------------------------------------------------------------------------------------------------
  gffread_convert_gff3_to_cds_fasta(file(params.new_assembly).getParent(),file(params.new_assembly).getName(),liftoff_annotations.out.liftoff_previous_annotations) // VALIDATED

  // -----------------------------------------------------------------------------------------------------------------------------------------------
  //         Run Salmon for strand inference and classify samples in three strand types : unstranded, stranded_forward and stranded_reverse
  // -----------------------------------------------------------------------------------------------------------------------------------------------
  salmon_index(gffread_convert_gff3_to_cds_fasta.out)

  salmon_strand_inference(trimming_fastq.out, salmon_index.out)

  def salmon_output_processed = salmon_strand_inference.out.map { sample_ID, library_layout, reads, strand_file ->
      def strand_info = file(strand_file).text.trim()
      return [sample_ID, library_layout, reads, strand_info]
  }

  salmon_output_processed.view { result ->
      println "DEBUG: salmon_output_processed -> ${result}"
      if (!result || result.isEmpty()) {
          println "ERROR: salmon_output_processed is empty!"
      } else {
          result.eachWithIndex { it, i -> println "DEBUG: Row[${i}] = ${it}" }
      }
  }

  def star_unstranded_out = Channel.empty()
  def hisat2_unstranded_out = Channel.empty()

  salmon_output_processed
      .collect()
      .view { list ->
          println "DEBUG: Checking if there are unstranded samples..."
          if (!list || list.isEmpty()) {
              println "WARNING: No samples found in `salmon_output_processed`."
          } else {
              println "DEBUG: Full list -> ${list}"
          }
      }
      .map { list -> list.any { it.size() > 3 && it[3] == "unstranded" } }
      .set { has_unstranded_samples }

  println "DEBUG: has_unstranded_samples = ${has_unstranded_samples}"

  // ----------------------------------------------------------------------------------------
  //                    Illumina short RNAseq reads alignment with STAR
  // ----------------------------------------------------------------------------------------
  star_genome_indices(file(params.new_assembly).getParent(),file(params.new_assembly).getName()) // VALIDATED

  // Always align stranded samples (stranded_forward and stranded_reverse)
  // Align ‘unstranded’ samples only if exist
  if (has_unstranded_samples) {
      star_alignment(star_genome_indices.out, salmon_output_processed) | collect | set { star_out }
  } else {
      println " WARNING : Skipping STAR alignment for unstranded samples, none found."
      star_alignment(star_genome_indices.out, salmon_output_processed
          .filter { sample_ID, library_layout, reads, strand_type -> strand_type in ["stranded_forward", "stranded_reverse"] }
      ) | collect | set { star_out }
  }

  star_out.view { result ->
    println "DEBUG: star_out -> ${result}"
  }

  // ----------------------------------------------------------------------------------------
  //                    Illumina short RNAseq reads alignment with HISAT2
  // ----------------------------------------------------------------------------------------
  hisat2_genome_indices(file(params.new_assembly).getParent(),file(params.new_assembly).getName())

  // Always align stranded samples (stranded_forward and stranded_reverse)
  // Align ‘unstranded’ samples only if exist
  if (has_unstranded_samples) {
      hisat2_alignment(hisat2_genome_indices.out, salmon_output_processed,
          file(params.new_assembly).getName()) | collect | set { hisat2_out }
  } else {
      println " WARNING : Skipping STAR alignment for unstranded samples, none found."
      hisat2_alignment(hisat2_genome_indices.out, salmon_output_processed
          .filter { sample_ID, library_layout, reads, strand_type -> strand_type in ["stranded_forward", "stranded_reverse"] },
          file(params.new_assembly).getName()) | collect | set { hisat2_out }
  }

  hisat2_out.view { result ->
      println "DEBUG: hisat2_out -> ${result}"
  }
 
  // ----------------------------------------------------------------------------------------
  //               Pacbio/Nanopore long RNAseq reads alignment with Minimap2  - OPTIONAL
  // ----------------------------------------------------------------------------------------

  // Check that samples_list_long_reads is empty or not before running minimap2-related processes
  minimap2_genome_indices(file(params.new_assembly).getParent(), file(params.new_assembly).getName())
  minimap2_alignment(minimap2_genome_indices.out, prepare_RNAseq_fastq_files_long.out) | collect

  minimap2_alignment
    .out
    .collect()
    .map { it[0] }
    .set { concat_minimap2_bams }

  has_long_reads.view { flag ->
      if (!flag) {
          println "No long reads detected, skipping Minimap2 genome indexing and alignment."
      } else {
          println "Long reads detected, running Minimap2 genome indexing and alignment."
      }
  }

  // ----------------------------------------------------------------------------------------
  //              transcriptome assembly with PsiCLASS on STAR alignments (short reads)
  // ----------------------------------------------------------------------------------------
/*  assembly_transcriptome_star_psiclass(star_stranded_out) | collect | set { star_psiclass_stranded_out }

  if (has_unstranded_samples) {
      assembly_transcriptome_star_psiclass(star_unstranded_out) | collect | set { star_psiclass_unstranded_out }
  } else {
      println "WARNING: No unstranded samples detected, skipping STAR/PsiCLASS for unstranded transcriptome."
  }*/

/*  // ----------------------------------------------------------------------------------------
  //        transcriptome assembly with Stringtie on STAR alignments (short reads)
  // ----------------------------------------------------------------------------------------
  // retrieve all the bam files to create a channel and launch StringTie one time per bam file
  // and then merge the transcriptomes

  star_alignment
  .out
  .collect()
  .flatten()
  .set{ concat_star_bams_stringtie } // VALIDATED
  assembly_transcriptome_star_stringtie(concat_star_bams_stringtie) | collect // VALIDATED
  assembly_transcriptome_star_stringtie
  .out
  .collect()
  .map { it[0] }
  .set{ concat_star_stringtie_annot } // VALIDATED
  Stringtie_merging_short_reads_STAR(concat_star_stringtie_annot) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //        transcriptome assembly with Stringtie on HISAT2 alignments (short reads)
  // ----------------------------------------------------------------------------------------
  // retrieve all the bam files to create a channel and launch StringTie one time per bam file
  // and then merge the transcriptomes

  hisat2_alignment
  .out
  .collect()
  .flatten()
  .set{ concat_hisat2_bams_stringtie } // VALIDATED
  assembly_transcriptome_hisat2_stringtie(concat_hisat2_bams_stringtie) | collect // VALIDATED
  assembly_transcriptome_hisat2_stringtie
  .out
  .collect()
  .map { it[0] }
  .set{ concat_hisat2_stringtie_annot } // VALIDATED
  Stringtie_merging_short_reads_hisat2(concat_hisat2_stringtie_annot) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //        gffcompare on HISAT2/Stringtie merged transcriptome assembly
  // ----------------------------------------------------------------------------------------

  gffcompare(Stringtie_merging_short_reads_hisat2.out)

  // ----------------------------------------------------------------------------------------
  //      transcriptome assembly with Stringtie on minimap2 alignments (long reads) - OPTIONAL
  // ----------------------------------------------------------------------------------------

  assembly_transcriptome_minimap2_stringtie(concat_minimap2_bams) | collect // VALIDATED
  assembly_transcriptome_minimap2_stringtie
    .out
    .collect()
    .map { it[0] }
    .set { concat_minimap2_stringtie_annot } // VALIDATED
  Stringtie_merging_long_reads(concat_minimap2_stringtie_annot)

  has_long_reads.view { flag ->
      if (!flag) {
          println "No long reads alignments detected, skipping transcriptome assembly and StringTie merging for long reads."
      } else {
          println "Long reads alignments detected, running transcriptome assembly and StringTie merging with long reads."
      }
  }

  // ----------------------------------------------------------------------------------------
  // -------------------------- Genome masking with EDTA ------------------------------------
  // ----------------------------------------------------------------------------------------

  EDTA(file(params.new_assembly).getParent(),file(params.new_assembly).getName()) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                                    BRAKER3 (AUGUSTUS/Genemark)
  // ----------------------------------------------------------------------------------------

  braker3_prediction(
      file(params.new_assembly).getParent(),
      file(params.new_assembly).getName(),
      file(params.protein_samplesheet).getParent(),
      file(params.protein_samplesheet).getName(),
      concat_star_bams_PsiCLASS
  )

  braker3_prediction_with_long_reads(
    file(params.new_assembly).getParent(),
    file(params.new_assembly).getName(),
    file(params.protein_samplesheet).getParent(),
    file(params.protein_samplesheet).getName(),
    concat_star_bams_PsiCLASS,
    concat_minimap2_bams
  )

  // ----------------------------------------------------------------------------------------
  //                                    Diamond2GO on proteins
  // ----------------------------------------------------------------------------------------
  // diamond2go(proteins_file)
*/
}
