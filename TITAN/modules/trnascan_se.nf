// tRNA gene prediction on the target genome with tRNAscan-SE.
process trnascan_se {
  label 'process_low'

  tag "tRNAscan-SE tRNA annotation"
  container params.container_trnascan
  publishDir "${params.output_dir}/additional_annotations/ncrna/trna", mode: 'copy'

  input:
    path(genome)
    path(trnascan_to_gff3)

  output:
    path "trnascan.out", emit: raw_table
    path "trnascan.struct", emit: structures
    path "trnascan.isotype", emit: isotypes
    path "trnascan.stats", emit: stats
    path "trna.gff3", emit: gff3
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    if [[ "${params.run_trnascan}" != "true" ]]; then
        printf "Sequence Name\\ttRNA #\\tBegin\\tEnd\\tType\\tAnticodon\\tIntron Begin\\tIntron End\\tScore\\n" > trnascan.out
        touch trnascan.struct trnascan.isotype trnascan.stats
        printf "##gff-version 3\\n" > trna.gff3
        printf '"%s":\n  trnascan_se: "skipped"\n  container: "%s"\n' \\
            "${task.process}" "${task.container}" > versions.yml
        exit 0
    fi

    if [[ ! -s "${genome}" ]]; then
        echo "ERROR: missing or empty target genome FASTA: ${genome}" >&2
        exit 1
    fi

    tRNAscan-SE -E \\
      -o trnascan.out \\
      -f trnascan.struct \\
      -s trnascan.isotype \\
      -m trnascan.stats \\
      --thread ${task.cpus} \\
      "${genome}"

    python3 "${trnascan_to_gff3}" trnascan.out > trna.gff3
    test -s trna.gff3

    version=\$(tRNAscan-SE -h 2>&1 | head -n 1 | sed 's/"/\\\\\\"/g')
    printf '"%s":\n  trnascan_se: "%s"\n  container: "%s"\n' \\
        "${task.process}" "\${version}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "Sequence Name\\ttRNA #\\tBegin\\tEnd\\tType\\tAnticodon\\tIntron Begin\\tIntron End\\tScore\\n" > trnascan.out
    printf "chrStub\\t1\\t10\\t80\\tAla\\tTGC\\t0\\t0\\t72.4\\n" >> trnascan.out
    printf ">chrStub.trna1 stub\\n" > trnascan.struct
    printf "Ala\\t1\\n" > trnascan.isotype
    printf "Total tRNAs\\t1\\n" > trnascan.stats
    printf "##gff-version 3\\nchrStub\\ttRNAscan-SE\\ttRNA\\t10\\t80\\t72.4\\t+\\t.\\tID=trnascan.chrStub.tRNA1;Name=chrStub.tRNA1;product=tRNA-Ala;anticodon=TGC\\n" > trna.gff3
    printf '"%s":\n  trnascan_se: "stub"\n' "${task.process}" > versions.yml
    """
}
