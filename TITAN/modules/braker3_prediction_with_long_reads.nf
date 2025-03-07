// AUGUSTUS can be run directly using the pipeline BRAKER3
process braker3_prediction_with_long_reads {

  tag "Executing BRAKER3/AUGUSTUS-Genemark prediction"
  container 'avelt/braker3:latest'
  containerOptions "--volume ${protein_samplesheet_path}:/protein_samplesheet_path --volume ${projectDir}/scripts:/scripts --volume ${projectDir}/work:/work --volume ${projectDir}/data/protein_data:/protein_path --volume ${genome_path}:/genome_path --volume ${params.output_dir}/intermediate_files/evidence_data/RNAseq_alignments/:/alignments --volume ${projectDir}:/outdir:z"
  cpus 4

  publishDir "${params.output_dir}", mode: 'copy'

  input:
    val(genome_path)
    val(genome)
    val(protein_samplesheet_path)
    val(protein_samplesheet_filename)
    val(bam_short)
    val(bam_long)

  output:
    path "augustus.hints.gff3", emit: augustus_gff
    path "genemark.gtf", emit: genemark_gtf
    path "genemark_supported.gtf", emit: genemark_supported_gtf
    path "braker.gff3", emit: braker_gff

  when:
  params.use_long_reads

  script:
    """
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running BRAKER3/AUGUSTUS-Genemark prediction"

    proteins=\$(/scripts/retrieve_proteins_for_braker.sh /protein_samplesheet_path/$protein_samplesheet_filename)
    bam_stranded_path=\$(/scripts/retrieve_path_bam_braker3.sh /alignments/STAR/stranded)
    bam_unstranded_path=\$(/scripts/retrieve_path_bam_braker3.sh /alignments/STAR/unstranded)

    if [ -z "\${bam_unstranded_path}" ]; then
        bam_merged="\${bam_stranded_path}"
    else
        bam_merged="\${bam_stranded_path},\${bam_unstranded_path}"
    fi

    bam_long_path=\$(/scripts/retrieve_path_bam_braker3.sh /alignments/minimap2)
    bam="\${bam_merged},\${bam_long_path}"

    CMD="/BRAKER-3.0.8/scripts/braker.pl --genome=/genome_path/$genome --bam=\${bam} \
    --prot_seq=\${proteins} \
    --threads=${task.cpus} --workingdir=\${PWD} --softmasking --gff3 \
    --PROTHINT_PATH=/ProtHint-2.6.0/bin/ --GENEMARK_PATH=/GeneMark-ETP --AUGUSTUS_CONFIG_PATH=/Augustus/config --TSEBRA_PATH=/TSEBRA/bin"
    echo "[\$DATE] Executing: \$CMD"

    /BRAKER-3.0.8/scripts/braker.pl --genome=/genome_path/$genome --bam=\${bam} \
    --prot_seq=\${proteins} \
    --threads=${task.cpus} --workingdir=\${PWD} --softmasking --gff3 \
    --PROTHINT_PATH=/ProtHint-2.6.0/bin/ --GENEMARK_PATH=/GeneMark-ETP --AUGUSTUS_CONFIG_PATH=/Augustus/config --TSEBRA_PATH=/TSEBRA/bin"
    # to test to decrease monoexon genes : --augustus_args="--singlestrand=true --alternatives-from-evidence=0"
    # --singlestrand=true: If enabled, only keeps predictions consistent with a single strand.
    # --alternatives-from-evidence=0: Disables alternative predictions based on hints, which can help filter out monoexons.
    cp Augustus/augustus.hints.gff3 .
    cp GeneMark-ETP/genemark.gtf .
    cp GeneMark-ETP/genemark_supported.gtf .
    """
}
