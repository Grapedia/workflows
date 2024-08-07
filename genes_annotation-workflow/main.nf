// Path to outdir
params.outdir = "${projectDir}/intermediate_files"

include { prepare_RNAseq_fastq_files } from "./modules/prepare_RNAseq_fastq_files"
include { trimming_fastq } from "./modules/trimming_fastq"
include { gmap_build_database } from "./modules/gmap_build_database"
include { gsnap_alignment } from "./modules/gsnap_alignment"
include { extracting_primary_mapping } from "./modules/extracting_primary_mapping"
// include { star_genome_indices } from "./modules/extracting_primary_mapping"
// include { star_alignment } from "./modules/extracting_primary_mapping"
include { assembly_transcriptome_stranded } from "./modules/assembly_transcriptome_stranded"
include { assembly_transcriptome_unstranded } from "./modules/assembly_transcriptome_unstranded"
include { conversion_gtf_gff3 } from "./modules/conversion_gtf_gff3"
include { split_proteins } from "./modules/split_proteins"
include { pblat_protein_alignment } from "./modules/pblat_protein_alignment"
include { exonerate_mapping } from "./modules/exonerate_mapping"
include { merge_exonerate_output } from "./modules/merge_exonerate_output"
include { filtering_and_conversion } from "./modules/filtering_and_conversion"
include { gtf_to_gff3 } from "./modules/gtf_to_gff3"
include { EDTA } from "./modules/EDTA"
include { liftoff_annotations } from "./modules/liftoff_annotations"
include { braker2_prediction_stranded } from "./modules/braker2_prediction_stranded"
include { braker2_prediction_unstranded } from "./modules/braker2_prediction_unstranded"
include { rename_braker2_gff_to_gff3 } from "./modules/rename_braker2_gff_to_gff3"
include { tRNAscan_SE } from "./modules/tRNAscan_SE"

Channel.fromPath( file(params.RNAseq_samplesheet) )
                    .splitCsv(header: true, sep: ',')
                    .map { row ->
                          def sample_ID = row.sampleID
                          def stranded_or_unstranded = row.stranded_or_unstranded
                          def SRA_or_FASTQ = row.SRA_or_FASTQ
                          def paired_or_single = row.paired_or_single

                          return [sample_ID, stranded_or_unstranded, SRA_or_FASTQ, paired_or_single]
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
  // -------------------------------- Evidence data analysis --------------------------------
  // ----------------------------------------------------------------------------------------

  // This step will generate for each assembly:
  //   * a transcriptome for stranded RNA-Seq data
  //   * a transcriptome for unstranded RNA-Seq data
  //   * an alignment file for protein sequences against the assembly
  // The transcriptomes will be generated based on RNA-Seq data and protein alignments ...
  // ... based on protein sequences

  // ----------------------------------------------------------------------------------------
  //                           Download/prepare RNAseq reads
  // ----------------------------------------------------------------------------------------
  prepare_RNAseq_fastq_files(samples_list) // VALIDATED
  trimming_fastq(prepare_RNAseq_fastq_files.out) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                           RNAseq reads alignment with GNSAP
  // ----------------------------------------------------------------------------------------
  gmap_build_database(params.assemblies_folder,params.new_assembly) // VALIDATED
  gsnap_alignment(gmap_build_database.out,trimming_fastq.out) // VALIDATED
  extracting_primary_mapping(gsnap_alignment.out) | collect // VALIDATED
  extracting_primary_mapping
  .out
  .collect()
  .map { it[0] }
  .set{ concat_in } // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                           RNAseq reads alignment with STAR
  // ----------------------------------------------------------------------------------------
  // star_genome_indices(params.assemblies_folder,params.new_assembly) // VALIDATED
  // star_alignment(star_genome_indices.out,trimming_fastq.out) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                           RNAseq reads alignment with HISAT2
  // ----------------------------------------------------------------------------------------


  // ----------------------------------------------------------------------------------------
  //              transcriptome assembly with PsiCLASS on GSNAP alignment
  // ----------------------------------------------------------------------------------------
  assembly_transcriptome_stranded(concat_in) // VALIDATED
  assembly_transcriptome_unstranded(concat_in) // VALIDATED
  conversion_gtf_gff3(params.assemblies_folder,params.new_assembly, assembly_transcriptome_stranded.out, assembly_transcriptome_unstranded.out) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                    Protein alignments with Pblat and Exonerate
  // ----------------------------------------------------------------------------------------
  split_proteins(protein_list) // VALIDATED
  pblat_protein_alignment(params.assemblies_folder,params.new_assembly, split_proteins.out.transpose()) // VALIDATED
  exonerate_mapping(params.assemblies_folder,params.new_assembly, pblat_protein_alignment.out)  | collect // VALIDATED
  exonerate_mapping
  .out
  .groupTuple()
  .set{ exonerate_map }
  merge_exonerate_output(exonerate_map) // VALIDATED
  filtering_and_conversion(merge_exonerate_output.out) // VALIDATED
  gtf_to_gff3(filtering_and_conversion.out) // VALIDATED

  // ----------------------------------------------------------------------------------------
  // -------------------------- Genome masking with EDTA ------------------------------------
  // ----------------------------------------------------------------------------------------

  EDTA(params.assemblies_folder,params.new_assembly) // VALIDATED

  // ----------------------------------------------------------------------------------------
  // --------------------------------- Ab initio predictions --------------------------------
  // ----------------------------------------------------------------------------------------

  // ----------------------------------------------------------------------------------------
  //                                Liftoff previous annotations
  // ----------------------------------------------------------------------------------------
  liftoff_annotations(params.assemblies_folder,params.new_assembly,params.annotations_folder,params.previous_assembly,params.previous_annotations) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                                    BRAKER2 (AUGUSTUS)
  // ----------------------------------------------------------------------------------------
  braker2_prediction_stranded(params.assemblies_folder,params.new_assembly,file(params.protein_samplesheet).getParent(),file(params.protein_samplesheet).getName(),concat_in)  // VALIDATED
  braker2_prediction_unstranded(params.assemblies_folder,params.new_assembly,file(params.protein_samplesheet).getParent(),file(params.protein_samplesheet).getName(),concat_in)  // VALIDATED
  // rename_braker2_gff_to_gff3(braker2_prediction_stranded.out,braker2_prediction_unstranded.out)

  // ----------------------------------------------------------------------------------------
  //                                  tRNAscan-SE annotation
  // ----------------------------------------------------------------------------------------

  tRNAscan_SE(params.assemblies_folder,params.new_assembly) // VALIDATED

}
