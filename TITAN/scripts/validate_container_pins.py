#!/usr/bin/env python3
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

REQUIRED_PARAMS = {
    "container_agat",
    "container_aegis",
    "container_braker3",
    "container_diamond2go",
    "container_eggnog_mapper",
    "container_edta",
    "container_egapx",
    "container_fastp",
    "container_gffcompare",
    "container_hisat2",
    "container_liftoff",
    "container_minimap2_samtools",
    "container_psiclass_samtools",
    "container_python",
    "container_salmon",
    "container_star",
    "container_stringtie",
}

PINNED_IMAGE = re.compile(r"^[\w./-]+@sha256:[0-9a-f]{64}$")


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def parse_nextflow_params():
    config = (ROOT / "nextflow.config").read_text()
    found = dict(re.findall(r"^\s*(container_[\w]+)\s*=\s*\"([^\"]+)\"", config, re.MULTILINE))
    missing = sorted(REQUIRED_PARAMS - set(found))
    if missing:
        return None, fail("missing container params: " + ", ".join(missing))

    invalid = [f"{name}={value}" for name, value in sorted(found.items()) if not PINNED_IMAGE.match(value)]
    if invalid:
        return None, fail("container params must use image@sha256 pins: " + ", ".join(invalid))

    for legacy_name, canonical_name in {"egapx_container": "container_egapx", "aegis_container": "container_aegis"}.items():
        if f"params.{legacy_name}" not in config or f"params.{canonical_name}" not in config:
            return None, fail(f"{legacy_name} compatibility alias for {canonical_name} is missing")
    return found, 0


def validate_modules():
    bad = []
    for path in sorted((ROOT / "modules").glob("*.nf")):
        text = path.read_text()
        for line_number, line in enumerate(text.splitlines(), 1):
            stripped = line.strip()
            if not stripped.startswith("container "):
                continue
            if stripped.startswith("container params.") or stripped.startswith('container "${params.'):
                continue
            bad.append(f"{path.relative_to(ROOT)}:{line_number}: {stripped}")
    if bad:
        return fail("module container directives must use central params:\n" + "\n".join(bad))
    return 0


def validate_dockerfiles():
    bad = []
    for path in sorted((ROOT / "dockerfiles").glob("**/Dockerfile")):
        for line_number, line in enumerate(path.read_text().splitlines(), 1):
            stripped = line.strip()
            if not stripped.startswith("FROM "):
                continue
            image = stripped.split()[1]
            if not PINNED_IMAGE.match(image):
                bad.append(f"{path.relative_to(ROOT)}:{line_number}: {stripped}")
    if bad:
        return fail("Dockerfile FROM directives must use image@sha256 pins:\n" + "\n".join(bad))
    return 0


def validate_latest_in_runtime_files():
    paths = [ROOT / "nextflow.config", *(ROOT / "modules").glob("*.nf")]
    offenders = []
    for path in paths:
        for line_number, line in enumerate(path.read_text().splitlines(), 1):
            if ":latest" in line:
                offenders.append(f"{path.relative_to(ROOT)}:{line_number}: {line.strip()}")
    if offenders:
        return fail("runtime files must not reference :latest images:\n" + "\n".join(offenders))
    return 0


def main():
    _, result = parse_nextflow_params()
    if result:
        return result
    for validator in (validate_modules, validate_dockerfiles, validate_latest_in_runtime_files):
        result = validator()
        if result:
            return result
    print("TITAN container pins OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
