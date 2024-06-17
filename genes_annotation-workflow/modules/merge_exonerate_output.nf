process merge_exonerate_output {

  tag "Merging exonerate output for $organism"
  container 'avelt/exonerate_bedtools_samtools:latest'
  containerOptions "--volume ${projectDir}/scripts/:/scripts --volume $params.outdir/evidence_data/protein_alignments_split/:/protein_alignments_split"
  publishDir "$params.outdir/evidence_data/protein_alignments_split/"
  cpus 4

  input:
    tuple val(organism), path(gff_file)

  output:
    tuple val(organism), file("*.gff")

  script:
    """
    if ls /protein_alignments_split/${organism}.gff 1> /dev/null 2>&1
    then
      cp /protein_alignments_split/${organism}.gff ${organism}.gff
    else
      /scripts/merge_exonerate_output.sh /protein_alignments_split/${organism} ${organism}
    fi
    """
}
