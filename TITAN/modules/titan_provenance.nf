process titan_provenance {
  label 'process_low'
  tag "TITAN evidence manifest"

  publishDir "${params.output_dir}/provenance", mode: 'copy'

  input:
    val(has_long_reads)
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

  output:
    path "evidence_manifest.json", emit: evidence_manifest
    path "versions.yml", emit: versions

  script:
    """
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

def records(label, value):
    paths = [item for item in str(value).split() if item]
    return [file_record(label, item) for item in paths]

manifest = {
    "schema_version": "titan.evidence_manifest.v1",
    "workflow": "TITAN",
    "has_long_reads": str("${has_long_reads}").lower() == "true",
    "params": {
        "output_dir": "${params.output_dir}",
        "egapx_version": "${params.egapx_version}",
        "egapx_revision": "${params.egapx_revision}",
        "egapx_container": "${params.container_egapx}",
        "egapx_data_version": "${params.egapx_data_version}",
        "aegis_version": "${params.aegis_version}",
        "aegis_container": "${params.container_aegis}",
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
    ] + records("long_reads_default", "${long_reads_default}") + records("long_reads_alt", "${long_reads_alt}"),
    "outputs": [
        file_record("aegis_gff", "${aegis_gff}"),
        file_record("aegis_proteins_all", "${aegis_proteins_all}"),
        file_record("aegis_proteins_main", "${aegis_proteins_main}"),
    ],
}

with open("evidence_manifest.json", "w", encoding="utf-8") as handle:
    json.dump(manifest, handle, indent=2, sort_keys=True)
    handle.write("\\n")

with open("versions.yml", "w", encoding="utf-8") as handle:
    handle.write('"TITAN:provenance":\\n')
    handle.write('  titan_manifest_schema: "titan.evidence_manifest.v1"\\n')
    handle.write('  egapx: "${params.egapx_version}"\\n')
    handle.write('  egapx_runner_revision: "${params.egapx_revision}"\\n')
    handle.write('  egapx_container: "${params.container_egapx}"\\n')
    handle.write('  aegis: "${params.aegis_version}"\\n')
    handle.write('  aegis_container: "${params.container_aegis}"\\n')
PY
    """
}
