#!/usr/bin/env python3
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = ROOT / "scripts" / "validate_final_annotation.py"


def write(path, text):
    path.write_text(text, encoding="utf-8")


def run(tmpdir, annotation_text):
    genome = tmpdir / "genome.fa"
    annotation = tmpdir / "annotation.gff3"
    proteins_all = tmpdir / "all.fa"
    proteins_main = tmpdir / "main.fa"
    write(genome, ">chr1\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n")
    write(annotation, annotation_text)
    write(proteins_all, ">prot1\nMKT\n")
    write(proteins_main, ">prot1\nMKT\n")
    return subprocess.run(
        [
            sys.executable,
            str(VALIDATOR),
            "--genome",
            str(genome),
            "--annotation",
            str(annotation),
            "--proteins-all",
            str(proteins_all),
            "--proteins-main",
            str(proteins_main),
            "--json-report",
            str(tmpdir / "report.json"),
            "--text-report",
            str(tmpdir / "report.txt"),
        ],
        text=True,
        capture_output=True,
    )


def expect_failure(name, tmpdir, annotation_text, expected):
    result = run(tmpdir, annotation_text)
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
    valid_gff = """##gff-version 3
chr1\tTITAN\tgene\t1\t39\t.\t+\t.\tID=gene1
chr1\tTITAN\tmRNA\t1\t39\t.\t+\t.\tID=gene1.t1;Parent=gene1
chr1\tTITAN\tCDS\t1\t39\t.\t+\t0\tID=gene1.t1.cds1;Parent=gene1.t1
"""
    failures = 0
    with tempfile.TemporaryDirectory() as raw:
        tmpdir = Path(raw)
        valid = run(tmpdir, valid_gff)
        if valid.returncode != 0:
            print("ERROR: valid final annotation fixture failed", file=sys.stderr)
            print(valid.stdout + valid.stderr, file=sys.stderr)
            return 1

        failures += expect_failure(
            "missing parent",
            tmpdir,
            "##gff-version 3\nchr1\tTITAN\tmRNA\t1\t39\t.\t+\t.\tID=tx1;Parent=missing_gene\n",
            "missing Parent target",
        )
        failures += expect_failure(
            "bad seqid",
            tmpdir,
            "##gff-version 3\nchrMissing\tTITAN\tgene\t1\t10\t.\t+\t.\tID=gene1\n",
            "is absent from genome FASTA",
        )
        failures += expect_failure(
            "bad CDS phase",
            tmpdir,
            "##gff-version 3\nchr1\tTITAN\tCDS\t1\t9\t.\t+\t.\tID=cds1\n",
            "CDS phase must be one of",
        )

    if failures:
        return 1
    print("TITAN final annotation validation tests OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
