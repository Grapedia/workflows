process filter_evidencemodeler_gff3 {

  tag "Executing EVM filtering on the following GFF3: $evidence_modeler_out"
  container 'avelt/agat_diamond:latest'
  containerOptions "--volume ${genome_path}:/genome_path --volume ${projectDir}/scripts/:/scripts --volume $NR_proteins_fasta_path:/NR_proteins_fasta_path --volume $uniprot_fasta_path:/uniprot_fasta_path"
  publishDir "$projectDir/FINAL_OUTPUT"
  cpus 4

  input:
    val(evidence_modeler_out)
    val(annotations_EVM_out)
    val(genome_path)
    val(genome)
    val(evm_at_least_2_ABINITIO_FINAL_gff3)
    val(evm_1_ABINITIO_FINAL_gff3)
    val(evm_evidencedata_only_FINAL_gff3)
    val(evm_1_ABINITIO_proteins_fasta)
    val(NR_proteins_fasta_path)
    val(NR_proteins_fasta_filename)
    val(uniprot_fasta_path)
    val(uniprot_fasta_filename)

  output:
    file("annotations.filtered.gff3")

  script:
    """
    /scripts/filter_annotations.sh -g /genome_path/$genome -i $evidence_modeler_out -r $annotations_EVM_out -o annotations.filtered.gff3 -a $evm_at_least_2_ABINITIO_FINAL_gff3 -b $evm_1_ABINITIO_FINAL_gff3 -c $evm_evidencedata_only_FINAL_gff3 -d $evm_1_ABINITIO_proteins_fasta -n /NR_proteins_fasta_path/$NR_proteins_fasta_filename -u /uniprot_fasta_path/$uniprot_fasta_filename -t ${task.cpus}
    """
}
