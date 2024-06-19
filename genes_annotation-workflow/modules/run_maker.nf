// For this step all the training and predictions with MAKER and SNAP are done using ...
// ... a bash script because MAKER doesn't allow to specify output path so it was easier ...
// ... to use bash
// MAKER will first use protein and RNA-Seq stranded transcripts to predict genes. SNAP ...
// ... is then trained and run a first time. With the output GFF file, the prediction tool ...
// ... is trained a second time because predicting again annotations

process run_maker {

  tag "Executing MAKER/SNAP prediction"
  container 'quay.io/biocontainers/maker:3.01.03--pl526hb8757ab_0'
  containerOptions "--volume ${params.outdir}:/outdir --volume ${projectDir}/scripts:/scripts --volume ${genome_path}:/genome_path --volume ${projectDir}/data/protein_data:/protein_path"
  publishDir "$params.outdir/MAKER"
  cpus 4

  input:
    val(genome_path)
    val(genome)
    val(RNAseq_stranded_transcriptome)
    val(protein_samplesheet)

  output:
    file("new_assembly_snap_second_prediction.all.gff")

  script:
    """
    proteins=\$(/scripts/retrieve_proteins_for_maker.sh $protein_samplesheet)
    /scripts/maker.sh -a /genome_path/${genome} \
    -t $RNAseq_stranded_transcriptome \
    -p \${proteins} \
    -o new_assembly -n ${task.cpus} -d /outdir
    """
}
