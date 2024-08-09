// Path to outdir
params.outdir = "${projectDir}/intermediate_files"

include { prepare_RNAseq_fastq_files } from "./modules/prepare_RNAseq_fastq_files"
include { trimming_fastq } from "./modules/trimming_fastq"
include { star_genome_indices } from "./modules/star_genome_indices"
include { star_alignment } from "./modules/star_alignment"
include { assembly_transcriptome_star_psiclass } from "./modules/assembly_transcriptome_star_psiclass"
include { conversion_gtf_gff3_star_psiclass } from "./modules/conversion_gtf_gff3_star_psiclass"
include { assembly_transcriptome_star_stringtie } from "./modules/assembly_transcriptome_star_stringtie"
include { EDTA } from "./modules/EDTA"
include { liftoff_annotations } from "./modules/liftoff_annotations"
include { braker3_prediction } from "./modules/braker3_prediction"
include { rename_braker3_gff_to_gff3 } from "./modules/rename_braker3_gff_to_gff3"
include { tRNAscan_SE } from "./modules/tRNAscan_SE"

Channel.fromPath( file(params.RNAseq_samplesheet) )
                    .splitCsv(header: true, sep: ',')
                    .map { row ->
                          def sample_ID = row.sampleID
                          def SRA_or_FASTQ = row.SRA_or_FASTQ
                          def paired_or_single = row.paired_or_single

                          return [sample_ID, SRA_or_FASTQ, paired_or_single]
                    }
                    .set{ samples_list }

Channel.fromPath( file(params.protein_samplesheet) )
                    .splitCsv(header: true, sep: ',')
                    .map { row ->
                          def organism = row.organism
                          def filename = row.filename
                          def braker3 = row.braker3

                          return [organism, filename, braker3]
                    }
                    .set{ protein_list }

workflow {

  // ----------------------------------------------------------------------------------------
  //                           Download/prepare RNAseq reads
  // ----------------------------------------------------------------------------------------
  prepare_RNAseq_fastq_files(samples_list) // VALIDATED
  trimming_fastq(prepare_RNAseq_fastq_files.out) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                           Illumina RNAseq reads alignment with STAR
  // ----------------------------------------------------------------------------------------
  star_genome_indices(file(params.new_assembly).getParent(),file(params.new_assembly).getName()) // VALIDATED
  star_alignment(star_genome_indices.out,trimming_fastq.out) | collect // VALIDATED
  // ----------------------------------------------------------------------------------------
  //                           Long RNAseq reads alignment with Minimap2
  // ----------------------------------------------------------------------------------------

  // ----------------------------------------------------------------------------------------
  //              transcriptome assembly with PsiCLASS on STAR alignments
  // ----------------------------------------------------------------------------------------
  // retrieve the first value to launch PsiClass assembly one time on all the bam files together
  star_alignment
  .out
  .collect()
  .map { it[0] }
  .set{ concat_star_bams_PsiCLASS } // VALIDATED
  assembly_transcriptome_star_psiclass(concat_star_bams_PsiCLASS) // VALIDATED
  conversion_gtf_gff3_star_psiclass(file(params.new_assembly).getParent(),file(params.new_assembly).getName(), assembly_transcriptome_star_psiclass.out) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //              transcriptome assembly with Stringtie on STAR alignments
  // ----------------------------------------------------------------------------------------
  // retrieve all the bam files to create a channel and launch StringTie one time per bam file
  star_alignment
  .out
  .collect()
  .flatten()
  .set{ concat_star_bams_stringtie } // VALIDATED
  concat_star_bams_stringtie.view()
  assembly_transcriptome_star_stringtie(concat_star_bams_stringtie) | collect
  // assembly_transcriptome_star_stringtie
  // .out
  // .collect()
  // .map { it[0] }
  // .set{ concat_stringtie_annot }
  // merge_annot_stringtie(concat_stringtie_annot)

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
  // braker3_prediction(file(params.new_assembly).getParent(),file(params.new_assembly).getName(),file(params.protein_samplesheet).getParent(),file(params.protein_samplesheet).getName(),concat_star_bams)  // VALIDATED
  // rename_braker3_gff_to_gff3(braker3_prediction.out)

  // ----------------------------------------------------------------------------------------
  //                                  tRNAscan-SE annotation
  // ----------------------------------------------------------------------------------------

  // tRNAscan_SE(file(params.new_assembly).getParent(),file(params.new_assembly).getName()) // VALIDATED

}
