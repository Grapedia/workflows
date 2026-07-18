process additional_annotations_provenance {
  label 'process_low'
  tag "TITAN additional annotations manifest"

  container params.container_python
  publishDir "${params.output_dir}/provenance", mode: 'copy'

  input:
    val(workflow_revision)
    val(workflow_commit_id)
    val(workflow_command_line)
    path(helixer_gff3)
    path(helixer_versions, stageAs: "module_versions/helixer_versions.yml")
    path(trnascan_gff3)
    path(trnascan_raw_table)
    path(trnascan_stats)
    path(trnascan_versions, stageAs: "module_versions/trnascan_versions.yml")
    path(rfam_gff3)
    path(rfam_tblout)
    path(rfam_search_log)
    path(rfam_versions, stageAs: "module_versions/rfam_versions.yml")

  output:
    path "additional_annotations_manifest.json", emit: manifest
    path "versions.yml", emit: versions

  script:
    """
set -euo pipefail
python3 - <<'PY'
import hashlib
import json
from pathlib import Path

def file_record(label, path):
    p = Path(path)
    if not p.exists():
        return {"label": label, "path": str(path), "present": False}
    if p.is_dir():
        return {"label": label, "path": str(p), "present": True, "type": "directory"}
    h = hashlib.sha256()
    with p.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return {
        "label": label,
        "path": str(p),
        "present": True,
        "type": "file",
        "size_bytes": p.stat().st_size,
        "sha256": h.hexdigest(),
    }

manifest = {
    "schema_version": "titan.additional_annotations_manifest.v1",
    "workflow": "TITAN",
    "workflow_metadata": {
        "revision": "${workflow_revision}",
        "commit_id": "${workflow_commit_id}",
        "command_line": "${workflow_command_line}",
    },
    "params": {
        "run_helixer": "${params.run_helixer}",
        "helixer_model": "${params.helixer_model}",
        "helixer_use_gpu": "${params.helixer_use_gpu}",
        "helixer_model_dir": "${params.helixer_model_dir}",
        "container_helixer": "${params.container_helixer}",
        "run_trnascan": "${params.run_trnascan}",
        "container_trnascan": "${params.container_trnascan}",
        "run_rfam": "${params.run_rfam}",
        "rfam_data_dir": "${params.rfam_data_dir}",
        "container_infernal": "${params.container_infernal}",
    },
    "outputs": [
        file_record("helixer_gff3", "${helixer_gff3}"),
        file_record("trnascan_gff3", "${trnascan_gff3}"),
        file_record("trnascan_raw_table", "${trnascan_raw_table}"),
        file_record("trnascan_stats", "${trnascan_stats}"),
        file_record("rfam_gff3", "${rfam_gff3}"),
        file_record("rfam_tblout", "${rfam_tblout}"),
        file_record("rfam_search_log", "${rfam_search_log}"),
    ],
    "module_versions": [
        file_record("helixer_versions", "${helixer_versions}"),
        file_record("trnascan_versions", "${trnascan_versions}"),
        file_record("rfam_versions", "${rfam_versions}"),
    ],
}

with open("additional_annotations_manifest.json", "w", encoding="utf-8") as handle:
    json.dump(manifest, handle, indent=2, sort_keys=True)
    handle.write("\\n")

with open("versions.yml", "w", encoding="utf-8") as handle:
    handle.write('"TITAN:additional_annotations_provenance":\\n')
    handle.write('  additional_annotations_manifest_schema: "titan.additional_annotations_manifest.v1"\\n')
    handle.write('  container: "${task.container}"\\n')
    handle.write('  helixer: "${params.run_helixer}"\\n')
    handle.write('  trnascan: "${params.run_trnascan}"\\n')
    handle.write('  rfam: "${params.run_rfam}"\\n')
PY
    """
}
