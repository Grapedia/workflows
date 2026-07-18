#!/usr/bin/env python3
"""Static quality checks for TITAN Nextflow modules."""

from __future__ import annotations

from pathlib import Path
import re
import sys
import json


ROOT = Path(__file__).resolve().parents[1]


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def module_blocks(text: str) -> tuple[str, str]:
    if "script:" not in text or "stub:" not in text:
        return "", ""
    script = text.split("script:", 1)[1].split("stub:", 1)[0]
    output = text.split("output:", 1)[1].split("script:", 1)[0] if "output:" in text else ""
    return output, script


main = (ROOT / "main.nf").read_text(encoding="utf-8")
if re.search(r"(?m)^\s*params\.[A-Za-z0-9_]+\s*=", main):
    fail("main.nf must not mutate params at runtime")

expected_main = """nextflow.enable.dsl = 2

include { TITAN } from './workflows/titan'

workflow {
    TITAN()
}
"""
if main != expected_main:
    fail("main.nf must stay minimal: DSL2 activation, TITAN include, workflow call")

for module_path in sorted((ROOT / "modules").glob("*.nf")):
    text = module_path.read_text(encoding="utf-8")
    if "\r" in text:
        fail(f"{module_path.relative_to(ROOT)} contains CRLF line endings")
    if re.search(r"(?m)[ \t]+$", text):
        fail(f"{module_path.relative_to(ROOT)} contains trailing whitespace")
    if re.search(r"emit\s+:", text):
        fail(f"{module_path.relative_to(ROOT)} uses spaced emit syntax")
    output, script = module_blocks(text)
    if re.search(r"(?m)^\s*file\s*\(", output):
        fail(f"{module_path.relative_to(ROOT)} uses file(...) in an output block")
    if script and "set -euo pipefail" not in script:
        fail(f"{module_path.relative_to(ROOT)} script block is missing set -euo pipefail")
    if output and "versions.yml" not in output:
        fail(f"{module_path.relative_to(ROOT)} output block is missing versions.yml")

for module_path in sorted((ROOT / "modules").glob("*.nf")):
    text = module_path.read_text(encoding="utf-8")
    forbidden_patterns = [
        'file("*.gtf")',
        'file("*.trimmed.fastq.gz")',
        'path("*-diamond*")',
        'path("${sample_ID}*")',
    ]
    for pattern in forbidden_patterns:
        if pattern in text:
            fail(f"{module_path.relative_to(ROOT)} contains broad output glob {pattern}")

for nextflow_path in [
    ROOT / "workflows" / "titan.nf",
    *sorted((ROOT / "subworkflows").glob("*.nf")),
    *sorted((ROOT / "modules").glob("*.nf")),
]:
    text = nextflow_path.read_text(encoding="utf-8")
    if "command.execute()" in text:
        fail(f"{nextflow_path.relative_to(ROOT)} runs local commands from workflow Groovy code")
    if ".ifEmpty([])" in text or "Channel.value([])" in text:
        fail(f"{nextflow_path.relative_to(ROOT)} uses an untyped empty list channel")
    if re.search(r"file\s*\([^)]*\)\.text", text):
        fail(f"{nextflow_path.relative_to(ROOT)} reads task files from workflow Groovy code")

for helper in ["run_aegis_merge.sh", "run_stringtie_transcriptome.sh", "clean_liftoff_gff3_for_agat.py"]:
    helper_path = ROOT / "scripts" / helper
    if not helper_path.exists():
        fail(f"shared helper script is missing: scripts/{helper}")
    if not helper_path.stat().st_mode & 0o111:
        fail(f"shared helper script is not executable: scripts/{helper}")

if not (ROOT / "docs" / "development" / "nextflow-dsl2-conventions.md").exists():
    fail("Nextflow DSL2 conventions document is missing")

generate_evidence = (ROOT / "subworkflows" / "generate_evidence_data.nf").read_text(encoding="utf-8")
if "params." in generate_evidence or "projectDir" in generate_evidence:
    fail("subworkflows/generate_evidence_data.nf must receive runtime files and options through take:")

schema_path = ROOT / "nextflow_schema.json"
if not schema_path.exists():
    fail("nextflow_schema.json is missing")

schema = json.loads(schema_path.read_text(encoding="utf-8"))
for required_param in [
    "output_dir",
    "new_assembly",
    "previous_assembly",
    "previous_annotations",
    "RNAseq_samplesheet",
    "RNAseq_data_dir",
    "protein_samplesheet",
    "egapx_paramfile",
    "container_egapx",
    "container_aegis",
    "container_eggnog_mapper",
    "container_helixer",
    "container_interproscan",
    "container_trnascan",
]:
    if required_param not in schema.get("properties", {}):
        fail(f"nextflow_schema.json is missing parameter: {required_param}")

print("Nextflow quality checks OK")
