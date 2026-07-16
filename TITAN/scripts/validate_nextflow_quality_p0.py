#!/usr/bin/env python3
"""Static P0 quality checks for TITAN Nextflow modules."""

from pathlib import Path
import re
import sys


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

for module_path in sorted((ROOT / "modules").glob("*.nf")):
    text = module_path.read_text(encoding="utf-8")
    output, script = module_blocks(text)
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
    if ".ifEmpty([])" in text or "Channel.value([])" in text:
        fail(f"{nextflow_path.relative_to(ROOT)} uses an untyped empty list channel")

print("P0 Nextflow quality checks OK")
