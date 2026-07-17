#!/usr/bin/env python3
import subprocess
import sys
import tempfile
import gzip
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
    "--run-eggnog-mapper",
    "false",
    "--eggnog-data-dir",
    str(ROOT / "data" / "eggnog"),
    "--run-helixer",
    "false",
    "--helixer-model-dir",
    str(ROOT / "data" / "helixer"),
    "--helixer-model",
    "land_plant",
    "--run-interproscan",
    "false",
    "--interproscan-data-dir",
    str(ROOT / "data" / "interproscan"),
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

    with tempfile.TemporaryDirectory(dir=INVALID) as temp_dir:
        temp_path = Path(temp_dir)
        fastq = "@sample.1 1 length=4\nACGT\n+sample.1 1 length=4\nIIII\n"
        with gzip.open(temp_path / "plus_header_1.fastq.gz", "wt") as handle:
            handle.write(fastq)
        with gzip.open(temp_path / "plus_header_2.fastq.gz", "wt") as handle:
            handle.write(fastq)
        samplesheet = temp_path / "rnaseq_plus_header.csv"
        samplesheet.write_text("sampleID,SRA_or_FASTQ,library_layout\nplus_header,FASTQ,paired\n")
        plus_header = run(
            replace_arg(
                replace_arg(BASE_ARGS, "--rnaseq-samplesheet", samplesheet),
                "--rnaseq-data-dir",
                temp_path,
            )
        )
        if plus_header.returncode != 0:
            print("ERROR: FASTQ plus-line headers should be accepted", file=sys.stderr)
            print(plus_header.stdout + plus_header.stderr, file=sys.stderr)
            failures += 1

    with tempfile.NamedTemporaryFile("w", suffix=".csv", dir=INVALID) as handle:
        handle.write("sampleID,SRA_or_FASTQ,library_layout\nnot_a_run,SRA,single\n")
        handle.flush()
        failures += expect_failure(
            "bad SRA accession",
            replace_arg(BASE_ARGS, "--rnaseq-samplesheet", handle.name),
            "SRA sampleID must be an SRR, ERR or DRR run accession",
        )

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
            "protein internal stop",
            replace_arg(BASE_ARGS, "--protein-samplesheet", INVALID / "protein_internal_stop.csv"),
            "internal stop codon",
        ),
        (
            "protein invalid character",
            replace_arg(BASE_ARGS, "--protein-samplesheet", INVALID / "protein_invalid_character.csv"),
            "invalid protein character",
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
        (
            "missing eggnog data dir",
            [
                *replace_arg(BASE_ARGS, "--run-eggnog-mapper", "true"),
                *["--eggnog-data-dir", str(INVALID / "missing_eggnog_dir")],
            ],
            "eggnog_data_dir does not exist or is not a directory",
        ),
        (
            "missing helixer model dir",
            [
                *replace_arg(BASE_ARGS, "--run-helixer", "true"),
                *["--helixer-model-dir", str(INVALID / "missing_helixer_dir")],
            ],
            "helixer_model_dir does not exist or is not a directory",
        ),
        (
            "missing helixer lineage model",
            [
                *replace_arg(BASE_ARGS, "--run-helixer", "true"),
                *["--helixer-model-dir", str(VALID)],
                *["--helixer-model", "fungi"],
            ],
            "no pre-fetched Helixer model for lineage 'fungi'",
        ),
        (
            "missing interproscan data dir",
            [
                *replace_arg(BASE_ARGS, "--run-interproscan", "true"),
                *["--interproscan-data-dir", str(INVALID / "missing_interproscan_dir")],
            ],
            "interproscan_data_dir does not exist or is not a directory",
        ),
        (
            "interproscan data dir missing member databases",
            [
                *replace_arg(BASE_ARGS, "--run-interproscan", "true"),
                *["--interproscan-data-dir", str(VALID)],
            ],
            "no pre-fetched InterProScan member database data",
        ),
    ]

    for name, args, expected in cases:
        failures += expect_failure(name, args, expected)

    with tempfile.NamedTemporaryFile("w", suffix=".csv", dir=INVALID) as handle:
        handle.write("organism,filename\nTerminalStop,test-data/minimal/invalid/protein_terminal_stop.fa\n")
        handle.flush()
        terminal_stop = run(replace_arg(BASE_ARGS, "--protein-samplesheet", handle.name))
        if terminal_stop.returncode != 0:
            print("ERROR: terminal protein stop should be accepted and cleaned before BRAKER3", file=sys.stderr)
            print(terminal_stop.stdout + terminal_stop.stderr, file=sys.stderr)
            failures += 1

    if failures:
        return 1
    print("TITAN input validation tests OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
