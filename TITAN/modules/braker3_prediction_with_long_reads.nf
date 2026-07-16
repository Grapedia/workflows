// AUGUSTUS can be run directly using the pipeline BRAKER3
process braker3_prediction_with_long_reads {
  label 'process_prediction'

  tag "Executing BRAKER3/AUGUSTUS-Genemark prediction"
  container params.container_braker3
  publishDir "${params.output_dir}", mode: 'copy'

  input:
    path(genome)
    path(protein_fastas)
    path(bam_short)
    path(bam_long)
    path(clean_protein_script)

  output:
    path "augustus.hints.gff3", emit: augustus_gff
    path "genemark.gtf", emit: genemark_gtf
    path "genemark_supported.gtf", emit: genemark_supported_gtf
    path "braker.gff3", emit: braker_gff
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    echo "[\$DATE] Running BRAKER3/AUGUSTUS-Genemark prediction"

    bam_merged=\$(printf '%s,' ${bam_short} | sed 's/,\$//')

    # Cleaned the protein fasta files for BRAKER3 -> simpler header and replace . and * by X
    cleaned_proteins=""
    for file in ${protein_fastas}; do
        cleaned="\$(basename "\${file%.*}").cleaned.fasta"
        python3 ${clean_protein_script} "\$file" "\$cleaned" "\${file%.*}"
        if [[ -z "\$cleaned_proteins" ]]; then
            cleaned_proteins="\$cleaned"
        else
            cleaned_proteins="\$cleaned_proteins,\$cleaned"
        fi
    done

    bam_long_path=\$(printf '%s,' ${bam_long} | sed 's/,\$//')
    bam="\${bam_merged},\${bam_long_path}"

    CMD="/BRAKER-3.0.8/scripts/braker.pl --genome=${genome} --bam=\${bam} \
    --prot_seq=\${cleaned_proteins} \
    --threads=${task.cpus} --workingdir=\${PWD} --softmasking --gff3 \
    --PROTHINT_PATH=/ProtHint-2.6.0/bin/ --GENEMARK_PATH=/GeneMark-ETP --AUGUSTUS_CONFIG_PATH=/Augustus/config --TSEBRA_PATH=/TSEBRA/bin"
    echo "[\$DATE] Executing: \$CMD"

    /BRAKER-3.0.8/scripts/braker.pl --genome=${genome} --bam=\${bam} \
    --prot_seq=\${cleaned_proteins} \
    --threads=${task.cpus} --workingdir=\${PWD} --softmasking --gff3 \
    --PROTHINT_PATH=/ProtHint-2.6.0/bin/ --GENEMARK_PATH=/GeneMark-ETP --AUGUSTUS_CONFIG_PATH=/Augustus/config --TSEBRA_PATH=/TSEBRA/bin
    # to test to decrease monoexon genes : --augustus_args="--singlestrand=true --alternatives-from-evidence=0"
    # --singlestrand=true: If enabled, only keeps predictions consistent with a single strand.
    # --alternatives-from-evidence=0: Disables alternative predictions based on hints, which can help filter out monoexons.
    cp Augustus/augustus.hints.gff3 .
    cp GeneMark-ETP/genemark.gtf .
    cp GeneMark-ETP/genemark_supported.gtf .
    test -s braker.gff3
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    printf "##gff-version 3\\nchr1\\tAUGUSTUS\\tgene\\t1\\t10\\t.\\t+\\t.\\tID=augustus_stub_gene\\n" > augustus.hints.gff3
    printf "chr1\\tGeneMark\\ttranscript\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"genemark_stub_gene\\"; transcript_id \\"genemark_stub_tx\\";\\n" > genemark.gtf
    cp genemark.gtf genemark_supported.gtf
    printf "##gff-version 3\\nchr1\\tBRAKER3\\tgene\\t1\\t10\\t.\\t+\\t.\\tID=braker_stub_gene\\n" > braker.gff3
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
