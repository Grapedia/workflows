process split_proteins {

  tag "Split proteins on ${organism}"
  container 'quay.io/biocontainers/biopython.convert:1.3.3--pyh5e36f6f_0'
  containerOptions "--volume ${projectDir}/data/protein_data:/protein_data --volume ${projectDir}/scripts/:/scripts"
  publishDir "$params.outdir/evidence_data/protein/fasta_splitted/${organism}"
  cpus 4

  input:
    tuple val(organism), val(filename), val(maker_braker2)

  output:
    tuple val(organism), file("*.fasta")

  script:
    """
    python3 /scripts/split_fasta.py /protein_data/$filename ${organism}_
    """
}
