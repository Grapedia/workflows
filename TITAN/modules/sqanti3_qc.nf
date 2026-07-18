process sqanti3_qc {
  label 'process_transcriptome'
  tag "SQANTI3 QC on ${source_label}"
  container params.container_sqanti3
  publishDir "${params.output_dir}/additional_annotations/sqanti3/${source_label}", mode: 'copy', saveAs: { filename ->
    if (filename.endsWith('_classification.txt') || filename.endsWith('_corrected.gtf') || filename.endsWith('_report.html') || filename.endsWith('_summary.tsv') || filename == 'versions.yml') {
      return filename
    }
    return null
  }

  input:
    val(source_label)
    path(isoforms_gtf)
    path(genome)
    path(reference_gff3)
    val(has_long_reads)

  output:
    path "${source_label}.sqanti3_classification.txt", emit: classification
    path "${source_label}.sqanti3_corrected.gtf", emit: corrected_gtf
    path "${source_label}.sqanti3_report.html", emit: html_report
    path "${source_label}.sqanti3_summary.tsv", emit: summary_tsv
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    classification="${source_label}.sqanti3_classification.txt"
    corrected_gtf="${source_label}.sqanti3_corrected.gtf"
    html_report="${source_label}.sqanti3_report.html"
    summary_tsv="${source_label}.sqanti3_summary.tsv"

    has_features() {
      awk 'NF && \$0 !~ /^#/ { found = 1; exit } END { exit found ? 0 : 1 }' "\$1"
    }

    write_empty_outputs() {
      local status="\$1"
      printf "isoform\\tstructural_category\\n" > "\${classification}"
      printf "# SQANTI3 %s for %s\\n" "\${status}" "${source_label}" > "\${corrected_gtf}"
      printf "<html><body>SQANTI3 %s for %s</body></html>\\n" "\${status}" "${source_label}" > "\${html_report}"
      printf "source\\tstatus\\ttotal_isoforms\\tFSM\\tISM\\tNIC\\tNNC\\tother\\tclassification\\tcorrected_gtf\\n" > "\${summary_tsv}"
      printf "%s\\t%s\\t0\\t0\\t0\\t0\\t0\\t0\\t%s\\t%s\\n" "${source_label}" "\${status}" "\${classification}" "\${corrected_gtf}" >> "\${summary_tsv}"
    }

    if [[ "${has_long_reads}" != "true" ]]; then
      write_empty_outputs "no_long_reads"
      printf '"%s":\n  sqanti3_qc: "no_long_reads"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
      exit 0
    fi

    if [[ "${params.run_sqanti3}" != "true" ]]; then
      write_empty_outputs "skipped"
      printf '"%s":\n  sqanti3_qc: "skipped"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
      exit 0
    fi

    if ! has_features "${isoforms_gtf}"; then
      write_empty_outputs "no_isoforms"
      printf '"%s":\n  sqanti3_qc: "no_isoforms"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
      exit 0
    fi

    # SQANTI3 requires GTF reference input.
    gffread "${reference_gff3}" -T -o reference.sqanti3.gtf

    # Provide libbz2.so.1 for gtfToGenePred when the image lacks that SONAME.
    sqanti3_libbz2_path="${params.sqanti3_libbz2_path}"
    if [[ -n "\${sqanti3_libbz2_path}" && "\${sqanti3_libbz2_path}" != "false" ]]; then
      if [[ ! -s "\${sqanti3_libbz2_path}" ]]; then
        echo "ERROR: SQANTI3 libbz2 compatibility file not found inside the container: \${sqanti3_libbz2_path}" >&2
        echo "Set --sqanti3_libbz2_path to the container-visible libbz2.so.1.* path, or false if the SQANTI3 image already provides libbz2.so.1." >&2
        exit 1
      fi
      mkdir -p .libfix
      ln -sf "\${sqanti3_libbz2_path}" .libfix/libbz2.so.1
      export LD_LIBRARY_PATH="\$(pwd)/.libfix:\${LD_LIBRARY_PATH:-}"
    fi

    sqanti3_qc.py \\
      --isoforms "${isoforms_gtf}" \\
      --refGTF reference.sqanti3.gtf \\
      --refFasta "${genome}" \\
      --dir . \\
      --output "${source_label}.sqanti3" \\
      --cpus ${task.cpus}

    if [[ ! -s "\${classification}" || ! -s "\${corrected_gtf}" ]]; then
      echo "ERROR: SQANTI3 finished without required outputs for ${source_label}." >&2
      echo "Expected non-empty files: \${classification}, \${corrected_gtf}" >&2
      exit 1
    fi
    [[ -s "\${html_report}" ]] || printf "<html><body>SQANTI3 report for %s</body></html>\\n" "${source_label}" > "\${html_report}"

    python3 - <<'PY'
from pathlib import Path

source = "${source_label}"
classification = Path("${source_label}.sqanti3_classification.txt")
summary = Path("${source_label}.sqanti3_summary.tsv")
counts = {"FSM": 0, "ISM": 0, "NIC": 0, "NNC": 0}
other = 0
total = 0
with classification.open("r", encoding="utf-8", errors="replace") as handle:
    header = handle.readline().rstrip("\\n").split("\\t")
    try:
        category_index = header.index("structural_category")
    except ValueError:
        category_index = None
    for line in handle:
        if not line.strip():
            continue
        total += 1
        fields = line.rstrip("\\n").split("\\t")
        category = fields[category_index] if category_index is not None and category_index < len(fields) else "other"
        if category in counts:
            counts[category] += 1
        else:
            other += 1
with summary.open("w", encoding="utf-8") as handle:
    handle.write("source\\tstatus\\ttotal_isoforms\\tFSM\\tISM\\tNIC\\tNNC\\tother\\tclassification\\tcorrected_gtf\\n")
    handle.write(f"{source}\\tcomplete\\t{total}\\t{counts['FSM']}\\t{counts['ISM']}\\t{counts['NIC']}\\t{counts['NNC']}\\t{other}\\t{classification}\\t{source}.sqanti3_corrected.gtf\\n")
PY

    sqanti3_qc.py --version 2>&1 | sed 's/^/  sqanti3: "/; s/\$/"/' | {
      printf '"%s":\\n' "${task.process}"
      cat
      printf '  container: "%s"\\n' "${task.container}"
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "isoform\\tstructural_category\\n%s_tx\\tFSM\\n" "${source_label}" > "${source_label}.sqanti3_classification.txt"
    printf "chrStub1\\tSQANTI3\\ttranscript\\t1\\t30\\t.\\t+\\t.\\tgene_id \\"%s_gene\\"; transcript_id \\"%s_tx\\";\\n" "${source_label}" "${source_label}" > "${source_label}.sqanti3_corrected.gtf"
    printf "<html><body>SQANTI3 stub for %s</body></html>\\n" "${source_label}" > "${source_label}.sqanti3_report.html"
    printf "source\\tstatus\\ttotal_isoforms\\tFSM\\tISM\\tNIC\\tNNC\\tother\\tclassification\\tcorrected_gtf\\n" > "${source_label}.sqanti3_summary.tsv"
    printf "%s\\tstub\\t1\\t1\\t0\\t0\\t0\\t0\\t%s.sqanti3_classification.txt\\t%s.sqanti3_corrected.gtf\\n" "${source_label}" "${source_label}" "${source_label}" >> "${source_label}.sqanti3_summary.tsv"
    printf '"%s":\n  sqanti3_qc: "stub"\n' "${task.process}" > versions.yml
    """
}

process sqanti3_qc_multiqc {
  label 'process_low'
  tag "SQANTI3 long-read isoform QC MultiQC summary"
  container params.container_python
  publishDir "${params.output_dir}/quality_report/sqanti3", mode: 'copy'

  input:
    path(stringtie_summary)
    path(flair_summary)

  output:
    path "sqanti3_long_read_isoform_qc_mqc.tsv", emit: multiqc_tsv
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    python3 - <<'PY'
from pathlib import Path

inputs = [Path("${stringtie_summary}"), Path("${flair_summary}")]
rows = []
for path in inputs:
    lines = [line.rstrip("\\n") for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(lines) < 2:
        continue
    header = lines[0].split("\\t")
    for line in lines[1:]:
        values = line.split("\\t")
        rows.append(dict(zip(header, values)))

with open("sqanti3_long_read_isoform_qc_mqc.tsv", "w", encoding="utf-8") as handle:
    handle.write("# id: titan_sqanti3_long_read_isoforms\\n")
    handle.write("# section_name: 'TITAN SQANTI3 long-read isoform QC'\\n")
    handle.write("# description: 'SQANTI3 structural-category summary for StringTie/Minimap2 long-read and FLAIR isoform assemblies.'\\n")
    handle.write("# plot_type: 'table'\\n")
    handle.write("Source\\tStatus\\tTotal isoforms\\tFSM\\tISM\\tNIC\\tNNC\\tOther\\n")
    for row in rows:
        handle.write("\\t".join([
            row.get("source", "unknown"),
            row.get("status", "unknown"),
            row.get("total_isoforms", "0"),
            row.get("FSM", "0"),
            row.get("ISM", "0"),
            row.get("NIC", "0"),
            row.get("NNC", "0"),
            row.get("other", "0"),
        ]) + "\\n")
PY
    printf '"%s":\n  sqanti3_qc_multiqc: "python-stdlib"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    cat > sqanti3_long_read_isoform_qc_mqc.tsv <<'EOF'
# id: titan_sqanti3_long_read_isoforms
# section_name: 'TITAN SQANTI3 long-read isoform QC'
# description: 'SQANTI3 structural categories for long-read isoforms.'
# plot_type: 'table'
Source	Status	Total isoforms	FSM	ISM	NIC	NNC	Other
stringtie_long_reads	stub	1	1	0	0	0	0
flair_isoforms	stub	1	1	0	0	0	0
EOF
    printf '"%s":\n  sqanti3_qc_multiqc: "stub"\n' "${task.process}" > versions.yml
    """
}
