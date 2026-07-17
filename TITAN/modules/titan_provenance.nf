process titan_provenance {
  label 'process_low'
  tag "TITAN evidence manifest"

  container params.container_python
  publishDir "${params.output_dir}/provenance", mode: 'copy'

  input:
    val(has_long_reads)
    val(output_dir)
    val(egapx_version)
    val(egapx_revision)
    val(egapx_container)
    val(egapx_data_version)
    val(aegis_version)
    val(aegis_container)
    val(workflow_revision)
    val(workflow_commit_id)
    val(workflow_nextflow_version)
    val(workflow_profile)
    val(workflow_command_line)
    path(new_assembly)
    path(previous_assembly)
    path(previous_annotations)
    path(rnaseq_samplesheet)
    path(protein_samplesheet)
    path(egapx_paramfile)
    path(masked_genome)
    path(liftoff_gff3)
    path(egapx_gff3)
    path(braker_augustus_gff)
    path(braker_genemark_gtf)
    path(star_stringtie_default_stranded)
    path(star_stringtie_alt_stranded)
    path(star_psiclass_stranded)
    path(star_psiclass_unstranded)
    path(star_stringtie_default_unstranded)
    path(star_stringtie_alt_unstranded)
    path(hisat2_stringtie_default_stranded)
    path(hisat2_stringtie_alt_stranded)
    path(hisat2_stringtie_default_unstranded)
    path(hisat2_stringtie_alt_unstranded)
    path(long_reads_default)
    path(long_reads_alt)
    path(aegis_gff)
    path(aegis_proteins_all)
    path(aegis_proteins_main)
    path(edta_versions, stageAs: "module_versions/edta_versions.yml")
    path(egapx_versions, stageAs: "module_versions/egapx_versions.yml")
    path(braker_versions, stageAs: "module_versions/braker_versions.yml")
    path(aegis_versions, stageAs: "module_versions/aegis_versions.yml")
    path(diamond2go_versions, stageAs: "module_versions/diamond2go_versions.yml")
    path(eggnog_versions, stageAs: "module_versions/eggnog_versions.yml")
    path(final_validation_versions, stageAs: "module_versions/final_validation_versions.yml")
    path(eggnog_annotations_all)
    path(eggnog_annotations_main)

  output:
    path "evidence_manifest.json", emit: evidence_manifest
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
    "schema_version": "titan.evidence_manifest.v1",
    "workflow": "TITAN",
    "has_long_reads": str("${has_long_reads}").lower() == "true",
    "workflow_metadata": {
        "revision": "${workflow_revision}",
        "commit_id": "${workflow_commit_id}",
        "nextflow_version": "${workflow_nextflow_version}",
        "profile": "${workflow_profile}",
        "command_line": "${workflow_command_line}",
    },
    "params": {
        "output_dir": "${output_dir}",
        "egapx_version": "${egapx_version}",
        "egapx_revision": "${egapx_revision}",
        "egapx_container": "${egapx_container}",
        "egapx_data_version": "${egapx_data_version}",
        "aegis_version": "${aegis_version}",
        "aegis_container": "${aegis_container}",
    },
    "inputs": [
        file_record("new_assembly", "${new_assembly}"),
        file_record("previous_assembly", "${previous_assembly}"),
        file_record("previous_annotations", "${previous_annotations}"),
        file_record("rnaseq_samplesheet", "${rnaseq_samplesheet}"),
        file_record("protein_samplesheet", "${protein_samplesheet}"),
        file_record("egapx_paramfile", "${egapx_paramfile}"),
    ],
    "evidence": [
        file_record("masked_genome", "${masked_genome}"),
        file_record("liftoff_gff3", "${liftoff_gff3}"),
        file_record("egapx_gff3", "${egapx_gff3}"),
        file_record("braker_augustus_gff", "${braker_augustus_gff}"),
        file_record("braker_genemark_gtf", "${braker_genemark_gtf}"),
        file_record("star_stringtie_default_stranded", "${star_stringtie_default_stranded}"),
        file_record("star_stringtie_alt_stranded", "${star_stringtie_alt_stranded}"),
        file_record("star_psiclass_stranded", "${star_psiclass_stranded}"),
        file_record("star_psiclass_unstranded", "${star_psiclass_unstranded}"),
        file_record("star_stringtie_default_unstranded", "${star_stringtie_default_unstranded}"),
        file_record("star_stringtie_alt_unstranded", "${star_stringtie_alt_unstranded}"),
        file_record("hisat2_stringtie_default_stranded", "${hisat2_stringtie_default_stranded}"),
        file_record("hisat2_stringtie_alt_stranded", "${hisat2_stringtie_alt_stranded}"),
        file_record("hisat2_stringtie_default_unstranded", "${hisat2_stringtie_default_unstranded}"),
        file_record("hisat2_stringtie_alt_unstranded", "${hisat2_stringtie_alt_unstranded}"),
        file_record("long_reads_default", "${long_reads_default}"),
        file_record("long_reads_alt", "${long_reads_alt}"),
    ],
    "outputs": [
        file_record("aegis_gff", "${aegis_gff}"),
        file_record("aegis_proteins_all", "${aegis_proteins_all}"),
        file_record("aegis_proteins_main", "${aegis_proteins_main}"),
        file_record("eggnog_annotations_all", "${eggnog_annotations_all}"),
        file_record("eggnog_annotations_main", "${eggnog_annotations_main}"),
    ],
    "module_versions": [
        file_record("edta_versions", "${edta_versions}"),
        file_record("egapx_versions", "${egapx_versions}"),
        file_record("braker_versions", "${braker_versions}"),
        file_record("aegis_versions", "${aegis_versions}"),
        file_record("diamond2go_versions", "${diamond2go_versions}"),
        file_record("eggnog_versions", "${eggnog_versions}"),
        file_record("final_validation_versions", "${final_validation_versions}"),
    ],
}

with open("evidence_manifest.json", "w", encoding="utf-8") as handle:
    json.dump(manifest, handle, indent=2, sort_keys=True)
    handle.write("\\n")

with open("versions.yml", "w", encoding="utf-8") as handle:
    handle.write('"TITAN:provenance":\\n')
    handle.write('  titan_manifest_schema: "titan.evidence_manifest.v1"\\n')
    handle.write('  container: "${task.container}"\\n')
    handle.write('  nextflow_version: "${workflow_nextflow_version}"\\n')
    handle.write('  workflow_revision: "${workflow_revision}"\\n')
    handle.write('  workflow_commit_id: "${workflow_commit_id}"\\n')
    handle.write('  egapx: "${egapx_version}"\\n')
    handle.write('  egapx_runner_revision: "${egapx_revision}"\\n')
    handle.write('  egapx_container: "${egapx_container}"\\n')
    handle.write('  aegis: "${aegis_version}"\\n')
    handle.write('  aegis_container: "${aegis_container}"\\n')
    handle.write('  eggnog_mapper: "${params.run_eggnog_mapper}"\\n')
    handle.write('  eggnog_data_dir: "${params.eggnog_data_dir}"\\n')
PY
    """
}
