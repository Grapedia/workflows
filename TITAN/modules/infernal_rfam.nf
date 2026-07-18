// Rfam covariance-model ncRNA annotation split by FASTA sequence.
process rfam_split_genome {
  label 'process_low'

  tag "Split target genome for parallel Infernal/Rfam"
  container params.container_python

  input:
    path(genome)

  output:
    path "rfam_genome_parts/*.fa", emit: fasta_parts
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    python3 - <<'PY'
from pathlib import Path
import re

outdir = Path("rfam_genome_parts")
outdir.mkdir(exist_ok=True)

def safe_name(raw):
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", raw.strip())
    return safe or "sequence"

current_id = None
handle = None
seen = {}
with open("${genome}", encoding="utf-8") as fasta:
    for line in fasta:
        if line.startswith(">"):
            if handle:
                handle.close()
            current_id = line[1:].strip().split()[0]
            base = safe_name(current_id)
            seen[base] = seen.get(base, 0) + 1
            suffix = "" if seen[base] == 1 else f".{seen[base]}"
            handle = open(outdir / f"{base}{suffix}.fa", "w", encoding="utf-8")
            handle.write(line)
        elif handle:
            handle.write(line)
if handle:
    handle.close()

if not any(outdir.glob("*.fa")):
    raise SystemExit("ERROR: no FASTA records found in ${genome}")
PY
    printf '"%s":\n  rfam_split_genome: "python-stdlib"\n  container: "%s"\n' \\
        "${task.process}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    mkdir -p rfam_genome_parts
    printf ">chrStub1\\nACGTACGTACGT\\n" > rfam_genome_parts/chrStub1.fa
    printf ">chrStub2\\nTGCATGCATGCA\\n" > rfam_genome_parts/chrStub2.fa
    printf '"%s":\n  rfam_split_genome: "stub"\n' "${task.process}" > versions.yml
    """
}

process infernal_rfam_search {
  label 'process_rfam'

  tag "Infernal/Rfam search on ${sequence_id}"
  container params.container_infernal

  input:
    tuple val(sequence_id), path(genome_part)

  output:
    tuple val(sequence_id), path("${sequence_id}.rfam_hits.tbl"), path("${sequence_id}.rfam_search.out"), emit: search_results
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    if [[ "${params.run_rfam}" != "true" ]]; then
        touch "${sequence_id}.rfam_hits.tbl" "${sequence_id}.rfam_search.out"
        printf '"%s":\n  infernal: "skipped"\n  rfam_data_dir: "%s"\n  container: "%s"\n' \\
            "${task.process}" "${params.rfam_data_dir}" "${task.container}" > versions.yml
        exit 0
    fi

    if [[ ! -s "${genome_part}" ]]; then
        echo "ERROR: missing or empty split genome FASTA: ${genome_part}" >&2
        exit 1
    fi

    if [[ -z "${params.rfam_data_dir}" || "${params.rfam_data_dir}" == "false" ]]; then
        echo "ERROR: rfam_data_dir must be provided when run_rfam=true" >&2
        exit 1
    fi

    if [[ ! -s "${params.rfam_data_dir}/Rfam.cm" || ! -s "${params.rfam_data_dir}/Rfam.clanin" ]]; then
        echo "ERROR: rfam_data_dir must contain Rfam.cm and Rfam.clanin" >&2
        exit 1
    fi

    cmsearch --cpu ${task.cpus} \\
      --tblout "${sequence_id}.rfam_hits.tbl" \\
      --cut_ga \\
      --rfam \\
      --nohmmonly \\
      --clanin "${params.rfam_data_dir}/Rfam.clanin" \\
      "${params.rfam_data_dir}/Rfam.cm" \\
      "${genome_part}" > "${sequence_id}.rfam_search.out"

    cmsearch_version=\$(cmsearch -h 2>&1 | head -n 2 | tr '\\n' ' ' | sed 's/"/\\\\\\"/g')
    printf '"%s":\n  infernal: "%s"\n  rfam_data_dir: "%s"\n  container: "%s"\n' \\
        "${task.process}" "\${cmsearch_version}" "${params.rfam_data_dir}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    cat > "${sequence_id}.rfam_hits.tbl" <<EOF
#target name        accession  query name  accession  mdl  mdl from  mdl to  seq from  seq to  strand  trunc  pass  gc  bias  score  E-value  inc  description of target
${sequence_id}      -          5S_rRNA     RF00001    cm   1         119     20        138     +       no     1     0.49 0.0   87.4   1e-20    !    stub 5S ribosomal RNA
EOF
    printf "# Infernal/Rfam stub for %s\\n" "${sequence_id}" > "${sequence_id}.rfam_search.out"
    printf '"%s":\n  infernal: "stub"\n' "${task.process}" > versions.yml
    """
}

process infernal_rfam_merge {
  label 'process_low'

  tag "Merge Infernal/Rfam split searches"
  container params.container_python
  publishDir "${params.output_dir}/additional_annotations/ncrna/rfam", mode: 'copy'

  input:
    path(rfam_search_results, stageAs: "rfam_search_results/*")
    path(rfam_tblout_to_gff3)

  output:
    path "rfam_hits.tbl", emit: tblout
    path "rfam_search.out", emit: search_log
    path "rfam_ncrna.gff3", emit: gff3
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail

    : > rfam_hits.tbl
    : > rfam_search.out

    while IFS= read -r -d '' tbl_file; do
        cat "\${tbl_file}" >> rfam_hits.tbl
    done < <(find rfam_search_results -name '*.rfam_hits.tbl' -print0 | sort -z)

    while IFS= read -r -d '' log_file; do
        printf "### %s ###\\n" "\${log_file}" >> rfam_search.out
        cat "\${log_file}" >> rfam_search.out
        printf "\\n" >> rfam_search.out
    done < <(find rfam_search_results -name '*.rfam_search.out' -print0 | sort -z)

    python3 "${rfam_tblout_to_gff3}" rfam_hits.tbl > rfam_ncrna.gff3
    test -s rfam_ncrna.gff3

    rfam_release="unknown"
    if [[ "${params.rfam_data_dir}" != "false" && -s "${params.rfam_data_dir}/Rfam.version" ]]; then
        rfam_release=\$(tr '\\n' ' ' < "${params.rfam_data_dir}/Rfam.version" | sed 's/"/\\\\\\"/g')
    fi
    printf '"%s":\n  infernal_rfam_merge: "python-stdlib"\n  rfam_data_dir: "%s"\n  rfam_release: "%s"\n  container: "%s"\n' \\
        "${task.process}" "${params.rfam_data_dir}" "\${rfam_release}" "${task.container}" > versions.yml
    """

  stub:
    """
    set -euo pipefail
    : > rfam_hits.tbl
    : > rfam_search.out
    while IFS= read -r -d '' tbl_file; do
        cat "\${tbl_file}" >> rfam_hits.tbl
    done < <(find rfam_search_results -name '*.rfam_hits.tbl' -print0 | sort -z)
    while IFS= read -r -d '' log_file; do
        printf "### %s ###\\n" "\${log_file}" >> rfam_search.out
        cat "\${log_file}" >> rfam_search.out
        printf "\\n" >> rfam_search.out
    done < <(find rfam_search_results -name '*.rfam_search.out' -print0 | sort -z)
    python3 "${rfam_tblout_to_gff3}" rfam_hits.tbl > rfam_ncrna.gff3
    printf '"%s":\n  infernal_rfam_merge: "stub"\n' "${task.process}" > versions.yml
    """
}
