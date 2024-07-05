// Path to outdir
params.outdir = "${projectDir}/intermediate_files"

include { prepare_RNAseq_fastq_files } from "./modules/prepare_RNAseq_fastq_files"
include { trimming_fastq } from "./modules/trimming_fastq"
include { gmap_build_database } from "./modules/gmap_build_database"
include { gsnap_alignment } from "./modules/gsnap_alignment"
include { extracting_primary_mapping } from "./modules/extracting_primary_mapping"
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
include { RepeatMasker } from "./modules/RepeatMasker"
include { liftoff_annotations } from "./modules/liftoff_annotations"
include { glimmerhmm_training } from "./modules/glimmerhmm_training"
include { split_fasta } from "./modules/split_fasta"
include { glimmerhmm_prediction } from "./modules/glimmerhmm_prediction"
// include { concat_glimmerhmm_prediction } from "./modules/concat_glimmerhmm_prediction"
// include { glimmerhmm_gff_to_gff3 } from "./modules/glimmerhmm_gff_to_gff3"
include { run_maker } from "./modules/run_maker"
// include { rename_maker_gff_to_gff3 } from "./modules/rename_maker_gff_to_gff3"
include { braker2_prediction_stranded } from "./modules/braker2_prediction_stranded"
include { braker2_prediction_unstranded } from "./modules/braker2_prediction_unstranded"
include { rename_braker2_gff_to_gff3 } from "./modules/rename_braker2_gff_to_gff3"
include { run_geneid } from "./modules/run_geneid"
// include { evidence_modeler } from "./modules/evidence_modeler"
// include { filter_evidencemodeler_gff3 } from "./modules/filter_evidencemodeler_gff3"
include { tRNAscan_SE } from "./modules/tRNAscan_SE"
// include { pasa } from "./modules/pasa"

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
                          def maker_braker2 = row.maker_braker2

                          return [organism, filename, maker_braker2]
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
  //                                Transcriptomes assembly
  // ----------------------------------------------------------------------------------------
  prepare_RNAseq_fastq_files(samples_list) // VALIDATED
  trimming_fastq(prepare_RNAseq_fastq_files.out) // VALIDATED
  gmap_build_database(params.assemblies_folder,params.new_assembly) // VALIDATED
  gsnap_alignment(gmap_build_database.out,trimming_fastq.out) // VALIDATED
  extracting_primary_mapping(gsnap_alignment.out) | collect // VALIDATED
  extracting_primary_mapping
  .out
  .collect()
  .map { it[0] }
  .set{ concat_in }
  assembly_transcriptome_stranded(concat_in) // VALIDATED
  assembly_transcriptome_unstranded(concat_in) // VALIDATED
  conversion_gtf_gff3(params.assemblies_folder,params.new_assembly, assembly_transcriptome_stranded.out, assembly_transcriptome_unstranded.out) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                                   Protein alignments
  // ----------------------------------------------------------------------------------------
  // split_proteins(protein_list) // VALIDATED
  // pblat_protein_alignment(params.assemblies_folder,params.new_assembly, split_proteins.out.transpose()) // VALIDATED
  // exonerate_mapping(params.assemblies_folder,params.new_assembly, pblat_protein_alignment.out)  | collect // VALIDATED
  // exonerate_mapping
  //  .out
  //  .groupTuple()
  //  .set{ exonerate_map }
  // merge_exonerate_output(exonerate_map) // VALIDATED
  // filtering_and_conversion(merge_exonerate_output.out) // VALIDATED
  // gtf_to_gff3(filtering_and_conversion.out) // VALIDATED

  // ----------------------------------------------------------------------------------------
  // -------------------------- Genome masking with EDTA/RepeatMasker -----------------------
  // ----------------------------------------------------------------------------------------

  // ----------------------------------------------------------------------------------------
  //                                        EDTA annotation
  // ----------------------------------------------------------------------------------------

  // EDTA(params.assemblies_folder,params.new_assembly) // ERROR : see https://github.com/oushujun/EDTA/issues/478
  // so we can't launch RepeatMasker, Maker and geneid for now

  // ----------------------------------------------------------------------------------------
  //                                        RepeatMasker
  // ----------------------------------------------------------------------------------------

  // RepeatMasker(params.assemblies_folder,params.new_assembly,EDTA.out.TElib_fasta)

  // ----------------------------------------------------------------------------------------
  // --------------------------------- Ab initio predictions --------------------------------
  // ----------------------------------------------------------------------------------------

  // 4 ab initio predictions will be done with different tools: SNAP (using MAKER), Augustus ...
  // ... (using BRAKER 2), GlimmerHMM and GeneID

  // ----------------------------------------------------------------------------------------
  //                                Liftoff previous annotations
  // ----------------------------------------------------------------------------------------
  liftoff_annotations(params.assemblies_folder,params.new_assembly,params.annotations_folder,params.previous_assembly,params.previous_annotations) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                                        GlimmerHMM
  // ----------------------------------------------------------------------------------------

  // GlimmerHMM is a gene finder based on Generalized Hidden Markov Model (GHMM) used to ...
  // ...predict annotations using JMM. This tool takes in input the RNA-Seq transcriptome ...
  // ... assembled with this pipeline and generates in output a TSV file containing the ...
  // ... exons location for each transcript (one exon per line and exons comming from ...
  // ... different transcripts are separated by a blank line)
  // glimmerhmm_training(params.assemblies_folder,params.new_assembly,conversion_gtf_gff3.out.stranded_gff3) // ERROR : ERROR 69: /usr/local/bin/score exited funny: 35584
  // split_fasta(params.assemblies_folder,params.new_assembly) // VALIDATED
  // glimmerhmm_prediction(split_fasta.out.flatten(),glimmerhmm_training.out) | collect
  // glimmerhmm_prediction
  // .out
  // .collect()
  // .map { it[0] }
  // .set{ glimmerhmm_pred }
  // concat_glimmerhmm_prediction(glimmerhmm_pred)
  // glimmerhmm_gff_to_gff3(concat_glimmerhmm_prediction.out)

  // ----------------------------------------------------------------------------------------
  //                                      MAKER (SNAP)
  // ----------------------------------------------------------------------------------------

  // MAKER is an annotation pipeline that integrates multiple predictor tools like SNAP
  // run_maker(RepeatMasker.out.masked_genome,conversion_gtf_gff3.out.stranded_gff3.parent,conversion_gtf_gff3.out.stranded_gff3.name,file(params.protein_samplesheet).getParent(),file(params.protein_samplesheet).getName())
  // rename_maker_gff_to_gff3(run_maker.out)

  // ----------------------------------------------------------------------------------------
  //                                    BRAKER2 (AUGUSTUS)
  // ----------------------------------------------------------------------------------------
  braker2_prediction_stranded(params.assemblies_folder,params.new_assembly,file(params.protein_samplesheet).getParent(),file(params.protein_samplesheet).getName(),concat_in)  // VALIDATED
  // braker2_prediction_unstranded(params.assemblies_folder,params.new_assembly,file(params.protein_samplesheet).getParent(),file(params.protein_samplesheet).getName(),concat_in)  // VALIDATED
  // rename_braker2_gff_to_gff3(braker2_prediction_stranded.out,braker2_prediction_unstranded.out)

  // ----------------------------------------------------------------------------------------
  //                                          GeneID
  // ----------------------------------------------------------------------------------------

  // Don' work on unmasked genome, so use the masked genome generated by RepeatMasker
  // run_geneid(RepeatMasker.out.masked_genome,file(params.geneid_param_file).getParent(),file(params.geneid_param_file).getName()) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                                        Evidence Modeler
  // ----------------------------------------------------------------------------------------

  // evidence_modeler(params.assemblies_folder,params.new_assembly,file(params.evm_config_file).getParent(),file(params.evm_config_file).getName(),run_geneid.out,glimmerhmm_gff_to_gff3.out,liftoff_annotations.out.liftoff_previous_annotations,rename_maker_gff_to_gff3.out,rename_braker2_gff_to_gff3.out.braker2_prediction_stranded,rename_braker2_gff_to_gff3.out.braker2_prediction_unstranded,conversion_gtf_gff3.out.stranded_gff3.parent,conversion_gtf_gff3.out.stranded_gff3.name,conversion_gtf_gff3.out.unstranded_gff3.parent,conversion_gtf_gff3.out.unstranded_gff3.name,EDTA.out.TE_annotations_gff3)

  // ----------------------------------------------------------------------------------------#
  //                              Evidence Modeler exonerate_filtering
  // ----------------------------------------------------------------------------------------#

  // filter_evidencemodeler_gff3(evidence_modeler.out.annotations_gff3,evidence_modeler.out.annotations_EVM_out,params.assemblies_folder,params.new_assembly,evidence_modeler.out.evm_at_least_2_ABINITIO_FINAL_gff3,evidence_modeler.out.evm_1_ABINITIO_FINAL_gff3,evidence_modeler.out.evm_evidencedata_only_FINAL_gff3,evidence_modeler.out.evm_1_ABINITIO_proteins_fasta,file(params.NR_proteins_fasta).getParent(),file(params.NR_proteins_fasta).getName(),file(params.uniprot_fasta).getParent(),file(params.uniprot_fasta).getName())

  // ----------------------------------------------------------------------------------------
  //                                  tRNAscan-SE annotation
  // ----------------------------------------------------------------------------------------

  tRNAscan_SE(params.assemblies_folder,params.new_assembly) // VALIDATED

  // ----------------------------------------------------------------------------------------
  //                                      PASA UTR annotation
  // ----------------------------------------------------------------------------------------

  // pasa(filter_evidencemodeler_gff3.out,params.assemblies_folder,params.new_assembly,conversion_gtf_gff3.out.stranded_fasta,file(params.pasa_config_file).getParent(),file(params.pasa_config_file).getName())

}
