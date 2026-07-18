process omark {
  label 'process_aegis'

  tag "Executing OMArk on final AEGIS main protein set"
  container params.container_omark
  publishDir "${params.output_dir}/quality_report/omark", mode: 'copy', saveAs: { filename ->
    if (filename in [
      'proteins_main.omamer',
      'proteins_main_detailed_summary.txt',
      'proteins_main_omark.sum',
      'omark_mqc.tsv',
      'versions.yml'
    ]) {
      return filename
    }
    return null
  }

  input:
    path(proteins_file_main)

  output:
    path "proteins_main.omamer", emit: omamer
    path "proteins_main_detailed_summary.txt", emit: detailed_summary
    path "proteins_main_omark.sum", emit: summary
    path "omark_mqc.tsv", emit: multiqc_tsv
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    write_empty_outputs() {
      local status="\$1"
      : > proteins_main.omamer
      printf "OMArk %s\\n" "\${status}" > proteins_main_detailed_summary.txt
      printf "OMArk %s\\n" "\${status}" > proteins_main_omark.sum
      cat > omark_mqc.tsv <<EOF
# id: titan_omark
# section_name: 'TITAN OMArk protein-set QC'
# description: 'OMArk completeness, consistency and contamination QC for the final AEGIS main protein set.'
# plot_type: 'table'
Metric	Value
Status	\${status}
Completeness	NA
Consistency	NA
Contamination	NA
EOF
    }

    if [[ "${params.run_omark}" != "true" ]]; then
      write_empty_outputs "skipped"
      printf '"%s":\n  omark: "skipped"\n  container: "%s"\n' "${task.process}" "${task.container}" > versions.yml
      exit 0
    fi

    if [[ ! -s "${proteins_file_main}" ]]; then
      echo "ERROR: missing or empty main proteins FASTA: ${proteins_file_main}" >&2
      exit 1
    fi

    if [[ -z "${params.omark_data_dir}" || "${params.omark_data_dir}" == "false" ]]; then
      echo "ERROR: omark_data_dir must be provided when run_omark=true" >&2
      exit 1
    fi

    omamer_db="${params.omark_data_dir}/omamer.h5"
    if [[ ! -s "\${omamer_db}" ]]; then
      echo "ERROR: OMArk/OMAmer database not found: \${omamer_db}" >&2
      exit 1
    fi

    omamer search \\
      --db "\${omamer_db}" \\
      --query "${proteins_file_main}" \\
      --out proteins_main.omamer \\
      --nthreads ${task.cpus}

    omark \\
      -f proteins_main.omamer \\
      -d "\${omamer_db}" \\
      -o omark_out

    cp omark_out/proteins_main_detailed_summary.txt proteins_main_detailed_summary.txt
    cp omark_out/proteins_main_omark.sum proteins_main_omark.sum

    test -s proteins_main_detailed_summary.txt
    test -s proteins_main_omark.sum

    python3 - <<'PY'
from pathlib import Path
import re

text = Path("proteins_main_omark.sum").read_text(encoding="utf-8", errors="replace")
patterns = {
    "Completeness": re.compile(r"complete(?:ness)?[^0-9]*([0-9]+(?:\\.[0-9]+)?%?)", re.IGNORECASE),
    "Consistency": re.compile(r"consistent|consistency", re.IGNORECASE),
    "Contamination": re.compile(r"contaminat[^0-9]*([0-9]+(?:\\.[0-9]+)?%?)", re.IGNORECASE),
}
values = {"Status": "complete", "Completeness": "see proteins_main_omark.sum", "Consistency": "see proteins_main_omark.sum", "Contamination": "see proteins_main_omark.sum"}
match = patterns["Completeness"].search(text)
if match:
    values["Completeness"] = match.group(1)
match = patterns["Contamination"].search(text)
if match:
    values["Contamination"] = match.group(1)
with open("omark_mqc.tsv", "w", encoding="utf-8") as handle:
    handle.write("# id: titan_omark\\n")
    handle.write("# section_name: 'TITAN OMArk protein-set QC'\\n")
    handle.write("# description: 'OMArk completeness, consistency and contamination QC for the final AEGIS main protein set.'\\n")
    handle.write("# plot_type: 'table'\\n")
    handle.write("Metric\\tValue\\n")
    for metric, value in values.items():
        handle.write(f"{metric}\\t{value}\\n")
PY

    {
      printf '"%s":\\n' "${task.process}"
      omark --version 2>&1 | sed 's/^/  omark: "/; s/\$/"/' || printf '  omark: "container-pinned"\\n'
      omamer --version 2>&1 | sed 's/^/  omamer: "/; s/\$/"/' || printf '  omamer: "container-pinned"\\n'
      printf '  omark_data_dir: "%s"\\n' "${params.omark_data_dir}"
      printf '  container: "%s"\\n' "${task.container}"
    } > versions.yml
    """

  stub:
    """
    set -euo pipefail
    printf "stub\\n" > proteins_main.omamer
    printf "OMArk stub detailed summary\\nCompleteness: 100%%\\nConsistency: pass\\nContamination: 0%%\\n" > proteins_main_detailed_summary.txt
    printf "Completeness: 100%%\\nConsistency: pass\\nContamination: 0%%\\n" > proteins_main_omark.sum
    cat > omark_mqc.tsv <<'EOF'
# id: titan_omark
# section_name: 'TITAN OMArk protein-set QC'
# description: 'OMArk completeness, consistency and contamination QC for the final AEGIS main protein set.'
# plot_type: 'table'
Metric	Value
Status	stub
Completeness	100%
Consistency	pass
Contamination	0%
EOF
    printf '"%s":\n  omark: "stub"\n  omamer: "stub"\n' "${task.process}" > versions.yml
    """
}
