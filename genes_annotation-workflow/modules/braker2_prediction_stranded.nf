// AUGUSTUS can be run directly using the pipeline BRAKER2
process braker2_prediction_stranded {

  tag "Executing BRAKER2/AUGUSTUS prediction on stranded data"
  container 'avelt/braker2_prothint_genemark:latest'
  containerOptions "--volume ${projectDir}/scripts:/scripts --volume ${genome_path}:/genome_path --volume $params.outdir/evidence_data/RNAseq_stranded/alignments/new_assembly:/alignments --volume ${projectDir}:/outdir"
  publishDir "$params.outdir/BRAKER2_RNAseq_stranded/"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    val(protein_samplesheet)
    val(sampleID)

  output:
    file("augustus.hints.gff3")

  script:
    """
    mkdir -p /outdir/TMP
    proteins=\$(/scripts/retrieve_proteins_for_maker.sh $protein_samplesheet)
    bam=\$(/scripts/retrieve_path_bam.sh /alignments)
    /BRAKER-2.1.6/scripts/braker.pl --genome=/genome_path/genome --bam=\${bam} \
    --prot_seq=\${proteins} \
    --cores=${task.cpus} --workingdir=/outdir/TMP --etpmode --softmasking --gff3 \
    --PROTHINT_PATH=/ProtHint-2.6.0/bin/ --GENEMARK_PATH=/GeneMark/
    """
}
