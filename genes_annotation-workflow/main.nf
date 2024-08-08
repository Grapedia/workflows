// Path to outdir
params.outdir = "${projectDir}/intermediate_files"

include { prepare_RNAseq_fastq_files } from "./modules/prepare_RNAseq_fastq_files"
include { trimming_fastq } from "./modules/trimming_fastq"
include { star_genome_indices } from "./modules/star_genome_indices"
include { star_alignment } from "./modules/star_alignment"
include { assembly_transcriptome_star_psiclass } from "./modules/assembly_transcriptome_star_psiclass"
include { conversion_gtf_gff3_star_psiclass } from "./modules/conversion_gtf_gff3_star_psiclass"
// include { split_proteins } from "./modules/split_proteins"
// include { pblat_protein_alignment } from "./modules/pblat_protein_alignment"
// include { exonerate_mapping } from "./modules/exonerate_mapping"
// include { merge_exonerate_output } from "./modules/merge_exonerate_output"
// include { filtering_and_conversion } from "./modules/filtering_and_conversion"
// include { gtf_to_gff3 } from "./modules/gtf_to_gff3"
include { EDTA } from "./modules/EDTA"
include { liftoff_annotations } from "./modules/liftoff_annotations"
include { braker2_prediction } from "./modules/braker2_prediction"
include { rename_braker2_gff_to_gff3 } from "./modules/rename_braker2_gff_to_gff3"
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
                          def braker2 = row.braker2

                          return [organism, filename, braker2]
                    }
                    .set{ protein_list }

workflow {

  // ----------------------------------------------------------------------------------------
  //                           Download/prepare RNAseq reads
  // ----------------------------------------------------------------------------------------
  prepare_RNAseq_fastq_files(samples_list) // VALIDATED
  trimming_fastq(prepare_RNAseq_fastq_files.out) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                           RNAseq reads alignment with STAR
  // ----------------------------------------------------------------------------------------
  star_genome_indices(file(params.new_assembly).getParent(),file(params.new_assembly).getName()) // VALIDATED
  star_alignment(star_genome_indices.out,trimming_fastq.out) | collect // VALIDATED
  star_alignment
  .out
  .collect()
  .map { it[0] }
  .set{ concat_star_bams } // VALIDATED

  // ----------------------------------------------------------------------------------------
  //              transcriptome assembly with PsiCLASS on STAR alignments
  // ----------------------------------------------------------------------------------------
  assembly_transcriptome_star_psiclass(concat_star_bams) // VALIDATED
  conversion_gtf_gff3_star_psiclass(file(params.new_assembly).getParent(),file(params.new_assembly).getName(), assembly_transcriptome_star_psiclass.out) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //              transcriptome assembly with Stringtie on STAR alignments
  // ----------------------------------------------------------------------------------------
  // assembly_transcriptome_star_stringtie(concat_star_bams) | collect
  // assembly_transcriptome_star_stringtie
  // .out
  // .collect()
  // .map { it[0] }
  // .set{ concat_stringtie_annot }
  // merge_annot_stringtie(concat_stringtie_annot)

  // ----------------------------------------------------------------------------------------
  //                    Protein alignments with Pblat and Exonerate
  // ----------------------------------------------------------------------------------------
  // split_proteins(protein_list) // VALIDATED
  // pblat_protein_alignment(file(params.new_assembly).getParent(),file(params.new_assembly).getName(), split_proteins.out.transpose()) // VALIDATED
  // exonerate_mapping(file(params.new_assembly).getParent(),file(params.new_assembly).getName(), pblat_protein_alignment.out)  | collect // VALIDATED
  // exonerate_mapping
  // .out
  // .groupTuple()
  // .set{ exonerate_map }
  // merge_exonerate_output(exonerate_map) // VALIDATED
  // filtering_and_conversion(merge_exonerate_output.out) // VALIDATED
  // gtf_to_gff3(filtering_and_conversion.out) // VALIDATED

  // ----------------------------------------------------------------------------------------
  // -------------------------- Genome masking with EDTA ------------------------------------
  // ----------------------------------------------------------------------------------------

  // EDTA(file(params.new_assembly).getParent(),file(params.new_assembly).getName()) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                                Liftoff previous annotations
  // ----------------------------------------------------------------------------------------
  // liftoff_annotations(file(params.new_assembly).getParent(),file(params.new_assembly).getName(),file(params.previous_assembly).getParent(),file(params.previous_assembly).getName(),file(params.previous_annotations).getParent(),file(params.previous_annotations).getName()) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                                    BRAKER2 (AUGUSTUS)
  // ----------------------------------------------------------------------------------------
  // braker2_prediction(file(params.new_assembly).getParent(),file(params.new_assembly).getName(),file(params.protein_samplesheet).getParent(),file(params.protein_samplesheet).getName(),concat_star_bams)  // VALIDATED
  // rename_braker2_gff_to_gff3(braker2_prediction.out)

  // ----------------------------------------------------------------------------------------
  //                                  tRNAscan-SE annotation
  // ----------------------------------------------------------------------------------------

  // tRNAscan_SE(file(params.new_assembly).getParent(),file(params.new_assembly).getName()) // VALIDATED

}
