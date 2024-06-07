// 5. Psiclass generates a GTF file, so we need to convert it in GFF3 format
process conversion_gtf_gff3 {

  tag "Convert GTF to GFF3 for transcriptome assembly"
  container 'quay.io/biocontainers/gffread:0.12.7--hdcf5f25_4'
  containerOptions "--volume $params.outdir/evidence_data/RNAseq_$stranded_or_unstranded/alignments/new_assembly:/alignments"
  publishDir "$params.outdir/evidence_data/transcriptomes/$stranded_or_unstranded"
  cpus 4

  input:

  output:
    file("$stranded_or_unstranded_vote.gtf")

  script:
    """
    ${projectDir}/scripts/psiclass_gff3_formatting.sh -g {input.assembly_file} -i {input.gtf} -o {output.gff3} -f {output.fasta} -s $stranded_or_unstranded
    """
}
