process filtering_and_conversion {

  tag "Executing filtering and GFF line extractions on the following protein sequences: ${organism}"
  container 'quay.io/biocontainers/biopython.convert:1.3.3--pyh5e36f6f_0'
  containerOptions "--volume ${projectDir}/scripts/:/scripts"
  publishDir "$params.outdir/evidence_data/protein_final_alignments"
  cpus 4

  input:
    tuple val(organism), path("gff_file")

  output:
    tuple val(organism), file("*.gff")

  script:
    """
    /scripts/filter_and_convert_exonerate_output.py -i $gff_file -o $organism.filtered.gff
    """
}
