process final_transcriptome_index {
  label 'process_index'
  tag "Final transcriptome Salmon index"
  container params.container_salmon
  publishDir "${params.output_dir}/quality_report/expression_validation", mode: 'copy', saveAs: { filename ->
    if (filename in ['final_transcripts.fasta', 'versions.yml']) {
      return filename
    }
    return null
  }

  input:
    path(genome)
    path(final_annotation_gff3)

  output:
    path "final_transcripts.fasta", emit: transcripts_fasta
    path "final_salmon_index", emit: index
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    if [[ "${params.run_expression_validation}" != "true" ]]; then
      : > final_transcripts.fasta
      mkdir -p final_salmon_index
      printf "expression validation skipped\\n" > final_salmon_index/versionInfo.json
      printf '"%s":\n  salmon: "skipped"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
      exit 0
    fi

    python3 - <<'PY'
from pathlib import Path

def parse_attrs(raw):
    attrs = {}
    for item in raw.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attrs[key] = value
    return attrs

def read_fasta(path):
    records = {}
    name = None
    seq = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(seq).upper()
            name = line[1:].split()[0]
            seq = []
        else:
            seq.append(line.strip())
    if name is not None:
        records[name] = "".join(seq).upper()
    return records

genome = read_fasta("${genome}")
tx_gene = {}
features = {}
fallback_features = {}
child_parent = {}
for line in Path("${final_annotation_gff3}").read_text().splitlines():
    if not line.strip() or line.startswith("#"):
        continue
    fields = line.split("\\t")
    if len(fields) != 9:
        continue
    seqid, source, feature_type, start, end, score, strand, phase, raw_attrs = fields
    attrs = parse_attrs(raw_attrs)
    feature_id = attrs.get("ID") or attrs.get("transcript_id")
    parent = attrs.get("Parent") or attrs.get("gene_id")
    if feature_id and parent:
        child_parent[feature_id] = parent.split(",")[0]
    if feature_type in {"mRNA", "transcript"}:
        transcript_id = attrs.get("ID") or attrs.get("transcript_id")
        gene_id = attrs.get("Parent") or attrs.get("gene_id")
        if transcript_id and gene_id:
            tx_gene[transcript_id] = gene_id.split(",")[0]
            fallback_features.setdefault(transcript_id, []).append((seqid, int(start), int(end), strand))
    if feature_type in {"exon", "CDS"} and parent:
        for transcript_id in parent.split(","):
            features.setdefault(transcript_id, []).append((seqid, int(start), int(end), strand))

for transcript_id, gene_id in list(tx_gene.items()):
    while gene_id in child_parent:
        gene_id = child_parent[gene_id]
    tx_gene[transcript_id] = gene_id

def revcomp(seq):
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]

with open("final_transcripts.fasta", "w", encoding="utf-8") as out:
    for transcript_id in sorted(tx_gene):
        parts = features.get(transcript_id) or fallback_features.get(transcript_id, [])
        if not parts:
            continue
        strand = parts[0][3]
        ordered = sorted(parts, key=lambda row: row[1], reverse=(strand == "-"))
        sequence = "".join(genome.get(seqid, "")[start - 1:end] for seqid, start, end, strand in ordered)
        if strand == "-":
            sequence = revcomp(sequence)
        if sequence:
            out.write(f">{transcript_id} gene_id={tx_gene[transcript_id]}\\n")
            out.write(sequence + "\\n")
PY

    if ! grep -q '^>' final_transcripts.fasta; then
      echo "ERROR: final transcriptome extraction produced no transcripts" >&2
      exit 1
    fi

    salmon index -t final_transcripts.fasta -i final_salmon_index -p ${task.cpus}
    salmon --version | sed 's/^/  salmon: "/; s/\$/"/' | {
      printf '"%s":\\n' "${task.process}"
      cat
      printf '  container: "%s"\\n' "${task.container}"
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf ">aegis_stub_tx gene_id=aegis_stub_gene\\nATGGCCATT\\n" > final_transcripts.fasta
    mkdir -p final_salmon_index
    printf "stub salmon index\\n" > final_salmon_index/versionInfo.json
    printf '"%s":\n  salmon: "stub"\n' "${task.process}" > versions.yml
    """
}

process final_expression_quant {
  label 'process_index'
  tag "Final transcriptome Salmon quant on ${sample_ID}"
  container params.container_salmon
  publishDir "${params.output_dir}/quality_report/expression_validation/quants", mode: 'copy', enabled: params.publish_intermediates, saveAs: { filename ->
    filename.endsWith('.trimmed.fastq.gz') ? null : filename
  }

  input:
    tuple val(sample_ID), val(library_layout), path(read_1), path(read_2)
    path(final_salmon_index)

  output:
    tuple val(sample_ID), path("${sample_ID}_quant"), emit: quant_dir
    path "${sample_ID}.expression_quant.log", emit: log
    path "${sample_ID}.expression_quant.versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    if [[ "${params.run_expression_validation}" != "true" ]]; then
      mkdir -p "${sample_ID}_quant"
      printf "Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\n" > "${sample_ID}_quant/quant.sf"
      printf "expression validation skipped\\n" > "${sample_ID}.expression_quant.log"
      printf '"%s":\n  salmon: "skipped"\n  container: "%s"\n' "${task.process}" "${task.container}" > "${sample_ID}.expression_quant.versions.yml"
      exit 0
    fi

    if [[ "${library_layout}" == "paired" ]]; then
      salmon quant -i "${final_salmon_index}" -l A -p ${task.cpus} \\
        -1 "${read_1}" -2 "${read_2}" -o "${sample_ID}_quant" --validateMappings \\
        2> "${sample_ID}.expression_quant.log"
    elif [[ "${library_layout}" == "single" ]]; then
      salmon quant -i "${final_salmon_index}" -l A -p ${task.cpus} \\
        -r "${read_1}" -o "${sample_ID}_quant" --validateMappings \\
        2> "${sample_ID}.expression_quant.log"
    else
      echo "ERROR: Unsupported library_layout '${library_layout}' for expression quantification" >&2
      exit 1
    fi

    test -s "${sample_ID}_quant/quant.sf"
    salmon --version | sed 's/^/  salmon: "/; s/\$/"/' | {
      printf '"%s":\\n' "${task.process}"
      cat
      printf '  container: "%s"\\n' "${task.container}"
    } > "${sample_ID}.expression_quant.versions.yml"
    """

  stub:
    """
    set -euo pipefail
    mkdir -p "${sample_ID}_quant"
    printf "Name\\tLength\\tEffectiveLength\\tTPM\\tNumReads\\n" > "${sample_ID}_quant/quant.sf"
    printf "aegis_stub_tx\\t9\\t9\\t10.0\\t10\\n" >> "${sample_ID}_quant/quant.sf"
    printf "stub expression quant\\n" > "${sample_ID}.expression_quant.log"
    printf '"%s":\n  salmon: "stub"\n' "${task.process}" > "${sample_ID}.expression_quant.versions.yml"
    """
}

process expression_support_summary {
  label 'process_low'
  tag "Final annotation expression support summary"
  container params.container_python
  publishDir "${params.output_dir}/quality_report/expression_validation", mode: 'copy'

  input:
    path(final_annotation_gff3)
    path(quant_dirs, stageAs: "quant_dirs/*")
    path(summarize_expression_support)

  output:
    path "expression_support_summary.json", emit: json_summary
    path "expression_support_summary_mqc.tsv", emit: multiqc_tsv
    path "gene_tpm_matrix.tsv", emit: tpm_matrix
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    if [[ "${params.run_expression_validation}" != "true" ]]; then
      cat > expression_support_summary.json <<EOF
{
  "min_tpm": ${params.expression_support_min_tpm},
  "samples": [],
  "total_genes": 0,
  "supported_genes": 0,
  "unsupported_genes": 0,
  "supported_gene_percent": 0.0,
  "unsupported_gene_percent": 0.0,
  "unsupported_gene_ids": [],
  "gene_tpm": {}
}
EOF
      cat > expression_support_summary_mqc.tsv <<EOF
# id: titan_expression_support
# section_name: 'TITAN expression support'
# description: 'Gene-level transcriptomic support from Salmon TPM on the final AEGIS transcriptome.'
# plot_type: 'table'
Metric	Value
Status	skipped
Total genes	0
Supported genes	0
Unsupported genes	0
Supported genes (%)	0.0
Minimum TPM	${params.expression_support_min_tpm}
EOF
      printf "gene_id\\n" > gene_tpm_matrix.tsv
      printf '"%s":\n  expression_support_summary: "skipped"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
      exit 0
    fi

    python3 "${summarize_expression_support}" \\
      --gff "${final_annotation_gff3}" \\
      --quant-dirs quant_dirs/* \\
      --min-tpm "${params.expression_support_min_tpm}" \\
      -o expression_support_summary.json \\
      --multiqc-tsv expression_support_summary_mqc.tsv \\
      --tpm-matrix gene_tpm_matrix.tsv

    printf '"%s":\n  expression_support_summary: "python-stdlib"\n  min_tpm: "%s"\n  container: "%s"\n' \\
      "${task.process}" "${params.expression_support_min_tpm}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    cat > expression_support_summary.json <<'EOF'
{
  "gene_tpm": {"aegis_stub_gene": {"synthetic_single": 10.0}},
  "min_tpm": 0.5,
  "samples": ["synthetic_single"],
  "supported_gene_percent": 100.0,
  "supported_genes": 1,
  "total_genes": 1,
  "unsupported_gene_ids": [],
  "unsupported_gene_percent": 0.0,
  "unsupported_genes": 0
}
EOF
    cat > expression_support_summary_mqc.tsv <<'EOF'
# id: titan_expression_support
# section_name: 'TITAN expression support'
# description: 'Gene-level transcriptomic support from Salmon TPM on the final AEGIS transcriptome.'
# plot_type: 'table'
Metric	Value
Total genes	1
Supported genes	1
Unsupported genes	0
Supported genes (%)	100.0
Minimum TPM	0.5
EOF
    printf "gene_id\\tsynthetic_single\\naegis_stub_gene\\t10.0\\n" > gene_tpm_matrix.tsv
    printf '"%s":\n  expression_support_summary: "stub"\n' "${task.process}" > versions.yml
    """
}
