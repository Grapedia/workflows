// AUGUSTUS can be run directly using the pipeline BRAKER3
process braker3_prediction {

  tag "Executing BRAKER3/AUGUSTUS-Genemark prediction"
  container 'avelt/braker3:latest'
  containerOptions "--volume ${protein_samplesheet_path}:/protein_samplesheet_path --volume ${projectDir}/scripts:/scripts --volume ${projectDir}/work:/work --volume ${projectDir}/data/protein_data:/protein_path --volume ${genome_path}:/genome_path --volume $params.outdir/evidence_data/RNAseq_alignments/STAR:/alignments --volume ${projectDir}:/outdir:z"
  publishDir "$params.outdir/BRAKER3/"
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
    proteins=\$(/scripts/retrieve_proteins_for_braker.sh /protein_samplesheet_path/$protein_samplesheet_filename)
    bam=\$(/scripts/retrieve_path_bam.sh /alignments)
    /BRAKER-3.0.8/scripts/braker.pl --genome=/genome_path/$genome --bam=\${bam} \
    --prot_seq=\${proteins} \
    --cores=${task.cpus} --workingdir=\${PWD} --etpmode --softmasking --gff3 \
    --PROTHINT_PATH=/ProtHint-2.6.0/bin/ --GENEMARK_PATH=/GeneMark/ --AUGUSTUS_CONFIG_PATH=/Augustus/config
    """
}
