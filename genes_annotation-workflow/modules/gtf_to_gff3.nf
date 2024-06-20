process gtf_to_gff3 {

  tag "Executing conversion gff to gff3 on the following protein sequences: $organism"
  container 'quay.io/biocontainers/evidencemodeler:1.1.1--0'
  publishDir "$params.outdir/evidence_data/protein_final_alignments/"
  cpus 4

  input:
    tuple val(organism), path(organism_filtered_gff_file)

  output:
    tuple val(organism), file("*.gff3")

  // replacement of /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/misc/Exonerate_to_evm_gff3.pl $organism_filtered_gff_file > ${organism}.gff3
  // by /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/misc/Exonerate_to_evm_gff3.pl $organism_filtered_gff_file > ${organism}.gff3
  script:
    """
    /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/misc/Exonerate_to_evm_gff3.pl $organism_filtered_gff_file > ${organism}.gff3
    """
}
