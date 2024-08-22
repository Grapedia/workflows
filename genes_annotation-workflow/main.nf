// Path to outdir
params.outdir = "${projectDir}/intermediate_files"

include { prepare_RNAseq_fastq_files_short } from "./modules/prepare_RNAseq_fastq_files_short"
include { prepare_RNAseq_fastq_files_long } from "./modules/prepare_RNAseq_fastq_files_long"
include { trimming_fastq } from "./modules/trimming_fastq"
include { star_genome_indices } from "./modules/star_genome_indices"
include { star_alignment } from "./modules/star_alignment"
include { minimap2_alignment } from "./modules/minimap2_alignment"
include { assembly_transcriptome_star_psiclass } from "./modules/assembly_transcriptome_star_psiclass"
include { assembly_transcriptome_star_stringtie } from "./modules/assembly_transcriptome_star_stringtie"
include { assembly_transcriptome_minimap2_stringtie } from "./modules/assembly_transcriptome_minimap2_stringtie"
include { Stringtie_merging_short_reads } from "./modules/Stringtie_merging_short_reads"
include { Stringtie_merging_long_reads } from "./modules/Stringtie_merging_long_reads"
include { EDTA } from "./modules/EDTA"
include { liftoff_annotations } from "./modules/liftoff_annotations"
include { braker3_prediction } from "./modules/braker3_prediction"
include { rename_braker3_gff_to_gff3 } from "./modules/rename_braker3_gff_to_gff3"
include { tRNAscan_SE } from "./modules/tRNAscan_SE"

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
  //                           Download/prepare RNAseq reads
  // ----------------------------------------------------------------------------------------
  prepare_RNAseq_fastq_files_short(samples_list_short_reads) // VALIDATED
  prepare_RNAseq_fastq_files_long(samples_list_long_reads) // VALIDATED
  // trimming with fastp in only done on Illumina short reads
  trimming_fastq(prepare_RNAseq_fastq_files_short.out) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                    Illumina short RNAseq reads alignment with STAR
  // ----------------------------------------------------------------------------------------
  star_genome_indices(file(params.new_assembly).getParent(),file(params.new_assembly).getName()) // VALIDATED
  star_alignment(star_genome_indices.out,trimming_fastq.out) | collect // VALIDATED
  // ----------------------------------------------------------------------------------------
  //               Pacbio/Nanopore long RNAseq reads alignment with Minimap2
  // ----------------------------------------------------------------------------------------

  minimap2_alignment(file(params.new_assembly).getParent(),file(params.new_assembly).getName(),prepare_RNAseq_fastq_files_long.out) | collect // VALIDATED
  minimap2_alignment
  .out
  .collect()
  .map { it[0] }
  .set{ concat_minimap2_bams } // VALIDATED

  // ----------------------------------------------------------------------------------------
  //              transcriptome assembly with PsiCLASS on STAR alignments (short reads)
  // ----------------------------------------------------------------------------------------
  // retrieve the first value to launch PsiClass assembly one time on all the bam files together

  star_alignment
  .out
  .collect()
  .map { it[0] }
  .set{ concat_star_bams_PsiCLASS } // VALIDATED
  assembly_transcriptome_star_psiclass(concat_star_bams_PsiCLASS) // VALIDATED

  // ----------------------------------------------------------------------------------------
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
  Stringtie_merging_short_reads(concat_star_stringtie_annot) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //      transcriptome assembly with Stringtie on minimap2 alignments (long reads)
  // ----------------------------------------------------------------------------------------

  assembly_transcriptome_minimap2_stringtie(minimap2_alignment.out) | collect // VALIDATED
  assembly_transcriptome_minimap2_stringtie
  .out
  .collect()
  .map { it[0] }
  .set{ concat_minimap2_stringtie_annot } // VALIDATED
  Stringtie_merging_long_reads(concat_minimap2_stringtie_annot) // VALIDATED

  // ----------------------------------------------------------------------------------------
  // -------------------------- Genome masking with EDTA ------------------------------------
  // ----------------------------------------------------------------------------------------

  // EDTA(file(params.new_assembly).getParent(),file(params.new_assembly).getName()) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                                Liftoff previous annotations
  // ----------------------------------------------------------------------------------------
  // liftoff_annotations(file(params.new_assembly).getParent(),file(params.new_assembly).getName(),file(params.previous_assembly).getParent(),file(params.previous_assembly).getName(),file(params.previous_annotations).getParent(),file(params.previous_annotations).getName()) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                                    BRAKER3 (AUGUSTUS/Genemark)
  // ----------------------------------------------------------------------------------------
  braker3_prediction(file(params.new_assembly).getParent(),file(params.new_assembly).getName(),file(params.protein_samplesheet).getParent(),file(params.protein_samplesheet).getName(),concat_star_bams_PsiCLASS,concat_minimap2_bams)
  // rename_braker3_gff_to_gff3(braker3_prediction.out)

  // ----------------------------------------------------------------------------------------
  //                                  tRNAscan-SE annotation
  // ----------------------------------------------------------------------------------------

  // tRNAscan_SE(file(params.new_assembly).getParent(),file(params.new_assembly).getName()) // VALIDATED

}
