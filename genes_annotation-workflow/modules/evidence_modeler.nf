// This step takes the following data to build consensus annotations based on weight ...
// ... given to each data:
//   * RNA-Seq stranded and unstranded transcriptomes
//   * protein alignments
//   * annotations predicted with the different ab initio prediction tools
//   * previous annotations lifted over the new assembly

process evidence_modeler {

  tag "Executing EvidenceModeler"
  container 'avelt/evidencemodeler_gffread:latest'
  containerOptions "--volume ${genome_path}:/genome_path --volume ${projectDir}/scripts/:/scripts --volume ${evm_config_file_path}:/evm_config_file_path --volume $params.outdir/evidence_data/protein_final_alignments/:/protein_final_alignments --volume $RNAseq_stranded_transcriptome_path:/RNAseq_stranded_transcriptome_path --volume $RNAseq_unstranded_transcriptome_path:/RNAseq_unstranded_transcriptome_path"
  publishDir "$projectDir/FINAL_OUTPUT"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    val(evm_config_file_path)
    val(evm_config_file_filename)
    val(run_geneid_out)
    val(glimmerhmm_predictions_out)
    val(liftoff_previous_annotations_out)
    val(maker_predictions_out)
    val(braker2_prediction_stranded)
    val(braker2_prediction_unstranded)
    val(RNAseq_stranded_transcriptome_path)
    val(RNAseq_stranded_transcriptome_filename)
    val(RNAseq_unstranded_transcriptome_path)
    val(RNAseq_unstranded_transcriptome_filename)

  output:
    path("annotations.gff3"), emit : annotations_gff3
    path("annotations.EVM.out"), emit : annotations_EVM_out
    path("evm.at_least_2_ABINITIO.FINAL.gff3"), emit : evm_at_least_2_ABINITIO_FINAL_gff3
    path("evm.1_ABINITIO.FINAL.gff3"), emit : evm_1_ABINITIO_FINAL_gff3
    path("evm.evidencedata_only.FINAL.gff3"), emit : evm_evidencedata_only_FINAL_gff3
    path("evm.1_ABINITIO.proteins.fasta"), emit : evm_1_ABINITIO_proteins_fasta

  script:
    """
    /scripts/evidencemodeler.sh -e \$PWD \
        -a /genome_path/$genome -w /evm_config_file_path/$evm_config_file_filename -i $run_geneid_out \
        -g $glimmerhmm_predictions_out -l $liftoff_previous_annotations_out -m $maker_predictions_out -b "$braker2_prediction_stranded braker2_prediction_unstranded" -r "" \
        -p /protein_final_alignments/*gff3 -t /RNAseq_unstranded_transcriptome_path/RNAseq_unstranded_transcriptome_filename -s /RNAseq_stranded_transcriptome_path/RNAseq_stranded_transcriptome_filename \
        -o annotations.gff3 -u annotations.EVM.out -c evm.at_least_2_ABINITIO.FINAL.gff3 -d evm.1_ABINITIO.FINAL.gff3 -f evm.evidencedata_only.FINAL.gff3 -j evm.1_ABINITIO.proteins.fasta
    """
}