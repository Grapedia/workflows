// AUGUSTUS can be run directly using the pipeline BRAKER2
process braker2_prediction_stranded {

  tag "Executing BRAKER2/AUGUSTUS prediction on stranded data"
  container 'avelt/braker2_prothint_genemark_augustus_bamtools_blast_samtools_diamond:latest'
  containerOptions "--volume ${protein_samplesheet_path}:/protein_samplesheet_path --volume ${projectDir}/scripts:/scripts --volume ${projectDir}/work:/work --volume ${projectDir}/data/protein_data:/protein_path --volume ${genome_path}:/genome_path --volume $params.outdir/evidence_data/RNAseq_stranded/alignments/new_assembly:/alignments --volume ${projectDir}:/outdir"
  publishDir "$params.outdir/BRAKER2_RNAseq_stranded/"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    val(protein_samplesheet_path)
    val(protein_samplesheet_filename)
    val(sampleID)

  output:
    file("augustus.hints.gff3")

  script:
    """
    mkdir -p /outdir/TMP
    chmod 777 /outdir/TMP
    proteins=\$(/scripts/retrieve_proteins_for_maker.sh /protein_samplesheet_path/$protein_samplesheet_filename)
    bam=\$(/scripts/retrieve_path_bam.sh /alignments)
    /BRAKER-2.1.6/scripts/braker.pl --genome=/genome_path/$genome --bam=\${bam} \
    --prot_seq=\${proteins} \
    --cores=${task.cpus} --workingdir=/outdir/TMP --etpmode --softmasking --gff3 \
    --PROTHINT_PATH=/ProtHint-2.6.0/bin/ --GENEMARK_PATH=/GeneMark/ --AUGUSTUS_CONFIG_PATH=/Augustus/config
    """
}
