process filtering_and_conversion {

  tag "Executing filtering and GFF line extractions on the following protein sequences: $organism"
  container 'avelt/exonerate_bedtools_samtools:latest'
  containerOptions "--memory=100g --volume ${projectDir}/scripts/:/scripts"
  publishDir "$params.outdir/evidence_data/protein_final_alignments/"
  cpus 4

  input:
    tuple val(organism), path(organism_gff_file)

  output:
    tuple val(organism), file("*.gff")

  script:
    """
    /scripts/filter_and_convert_exonerate_output.py -i $organism_gff_file -o ${organism}.filtered.gff
    """
}
