process agat_convert_gff3_to_cds_fasta {

  tag "Executing agat to convert liftoff gff3 to cds.fasta"
  container 'quay.io/biocontainers/agat:1.2.0--pl5321hdfd78af_0'
  publishDir "${params.output_dir}/intermediate_files/liftoff/gff3_to_cds_fasta"
  cpus 4

  input:
    path(genome)
    path(liftoff_gff3)

  output:
    file("${genome.simpleName}.CDS.fasta.gz")

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running agat_convert_gff3_to_cds_fasta on ${genome}"

    # formatting gff3 for AGAT (mandatory for PN40024 V4.3 gff3 version from grapedia)
    grep "obsolete=true" ${liftoff_gff3} | grep "Deleted" | awk -F'\t' '{print \$9}' | sed 's/.*ID=\\([^;]*\\);.*/\\1/' > to_remove.txt || true
    grep -F -f to_remove.txt ${liftoff_gff3} | grep "mRNA" | awk -F'\t' '{print \$9}' | sed 's/.*Parent=\\([^;]*\\);.*/\\1/' >> to_remove.txt || true
    grep -v -F -f to_remove.txt ${liftoff_gff3} > cleaned.gff3
    awk 'BEGIN{FS=OFS="\t"} /^#/ || NF==9' cleaned.gff3 > cleaned.OK.gff3
    fold -w 80 ${genome} > reformatted.fa
    CMD="agat_sp_extract_sequences.pl -g cleaned.OK.gff3 -f reformatted.fa -o ${genome.simpleName}.CDS.fasta"
    echo "[\$DATE] Executing: \$CMD"
    agat_sp_extract_sequences.pl -g cleaned.OK.gff3 -f reformatted.fa -o ${genome.simpleName}.CDS.fasta
    gzip ${genome.simpleName}.CDS.fasta
    """

  stub:
    """
    printf ">stub_cds\\nATGGCC\\n" | gzip -c > ${genome.simpleName}.CDS.fasta.gz
    """
}
