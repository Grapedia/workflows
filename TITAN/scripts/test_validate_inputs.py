#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VALID = ROOT / "test-data" / "minimal" / "valid"
INVALID = ROOT / "test-data" / "minimal" / "invalid"
VALIDATOR = ROOT / "scripts" / "validate_inputs.py"


BASE_ARGS = [
    sys.executable,
    str(VALIDATOR),
    "--project-dir",
    str(ROOT),
    "--new-assembly",
    str(VALID / "target.fa"),
    "--previous-assembly",
    str(VALID / "reference.fa"),
    "--previous-annotations",
    str(VALID / "reference.gff3"),
    "--rnaseq-samplesheet",
    str(VALID / "rnaseq_samplesheet.csv"),
    "--rnaseq-data-dir",
    str(VALID / "rnaseq"),
    "--protein-samplesheet",
    str(VALID / "protein_samplesheet.csv"),
    "--egapx-paramfile",
    str(VALID / "input_egapx.yaml"),
    "--egapx-executor",
    "docker",
    "--psiclass-vd",
    "5.0",
    "--psiclass-c",
    "0.03",
]


def run(args):
    return subprocess.run(args, cwd=ROOT, text=True, capture_output=True)


def replace_arg(args, flag, value):
    updated = list(args)
    updated[updated.index(flag) + 1] = str(value)
    return updated


def expect_failure(name, args, expected):
    result = run(args)
    output = result.stdout + result.stderr
    if result.returncode == 0:
        print(f"ERROR: {name}: expected failure but command succeeded", file=sys.stderr)
        return 1
    if expected not in output:
        print(f"ERROR: {name}: expected message fragment not found: {expected}", file=sys.stderr)
        print(output, file=sys.stderr)
        return 1
    return 0


def main():
    failures = 0
    valid = run(BASE_ARGS)
    if valid.returncode != 0:
        print("ERROR: valid fixture command failed", file=sys.stderr)
        print(valid.stdout + valid.stderr, file=sys.stderr)
        return 1

    cases = [
        (
            "bad RNA-seq layout",
            replace_arg(BASE_ARGS, "--rnaseq-samplesheet", INVALID / "rnaseq_bad_layout.csv"),
            "library_layout must be one of",
        ),
        (
            "missing RNA-seq FASTQ",
            replace_arg(BASE_ARGS, "--rnaseq-samplesheet", INVALID / "rnaseq_missing_fastq.csv"),
            "RNA-seq FASTQ for sample missing_fastq does not exist",
        ),
        (
            "duplicate RNA-seq sample",
            replace_arg(BASE_ARGS, "--rnaseq-samplesheet", INVALID / "rnaseq_duplicate_sample.csv"),
            "duplicate sampleID",
        ),
        (
            "missing protein FASTA",
            replace_arg(BASE_ARGS, "--protein-samplesheet", INVALID / "protein_missing_file.csv"),
            "protein FASTA for MissingProtein does not exist",
        ),
        (
            "bad previous GFF3 coordinates",
            replace_arg(BASE_ARGS, "--previous-annotations", INVALID / "bad_coordinates.gff3"),
            "invalid coordinates",
        ),
        (
            "missing EGAPx genome",
            replace_arg(BASE_ARGS, "--egapx-paramfile", INVALID / "egapx_missing_genome.yaml"),
            "EGAPx YAML must define 'genome'",
        ),
        (
            "bad EGAPx executor",
            replace_arg(BASE_ARGS, "--egapx-executor", "bad_executor"),
            "egapx_executor must be one of",
        ),
    ]

    for name, args, expected in cases:
        failures += expect_failure(name, args, expected)

    if failures:
        return 1
    print("TITAN input validation tests OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
