process gtf_to_gff3 {

  tag "Executing conversion gff to gff3 on the following protein sequences: ${organism}"
  container 'quay.io/biocontainers/evidencemodeler:1.1.1--hdfd78af_3'
  publishDir "$params.outdir/evidence_data/protein_final_alignments"
  cpus 4

  input:
    tuple val(organism), path("gff_file")

  output:
    tuple val(organism), file("*.gff")

  script:
    """
    /usr/local/opt/evidencemodeler-1.1.1/EvmUtils/misc/exonerate_gff_to_alignment_gff3.pl $gff_file prot > ${organism}.gff3
    """
}
