#!/usr/bin/env python3
import csv
import gzip
import hashlib
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
DATASET = ROOT / "test-data" / "minimal"
VALID = DATASET / "valid"
CHECKSUMS = DATASET / "checksums.sha256"


REQUIRED_VALID_FILES = [
    "reference.fa",
    "target.fa",
    "reference.gff3",
    "input_egapx.yaml",
    "rnaseq_samplesheet.csv",
    "protein_samplesheet.csv",
    "rnaseq/synthetic_single.fastq.gz",
    "rnaseq/synthetic_paired_1.fastq.gz",
    "rnaseq/synthetic_paired_2.fastq.gz",
    "rnaseq/synthetic_long.fastq.gz",
    "proteins/synthetic_core.fa",
    "proteins/synthetic_related.fa",
    "evidence/assembly_masked.EDTA.fasta",
    "evidence/liftoff_previous_annotations.gff3",
    "evidence/augustus.hints.gff3",
    "evidence/genemark.gtf",
    "evidence/merged_star_stringtie_stranded_default.gtf",
    "evidence/merged_star_stringtie_stranded_alt.gtf",
    "evidence/merged_star_stringtie_unstranded_default.gtf",
    "evidence/merged_star_stringtie_unstranded_alt.gtf",
    "evidence/merged_star_psiclass_stranded.gtf",
    "evidence/merged_star_psiclass_unstranded.gtf",
    "evidence/merged_minimap2_stringtie_long_reads_default.gtf",
    "evidence/merged_minimap2_stringtie_long_reads_alt.gtf",
]


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def read_fasta(path):
    seqs = {}
    current = None
    chunks = []
    with path.open() as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    seqs[current] = "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current:
        seqs[current] = "".join(chunks)
    return seqs


def validate_gff3(path, seqids):
    ids = set()
    parents = []
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                return fail(f"{path}: line {line_number} has {len(fields)} columns, expected 9")
            seqid, _, _, start, end, _, _, _, attributes = fields
            if seqid not in seqids:
                return fail(f"{path}: line {line_number} references unknown seqid {seqid}")
            if int(start) < 1 or int(end) < int(start) or int(end) > len(seqids[seqid]):
                return fail(f"{path}: line {line_number} has invalid coordinates")
            attrs = dict(item.split("=", 1) for item in attributes.split(";") if "=" in item)
            if "ID" in attrs:
                ids.add(attrs["ID"])
            if "Parent" in attrs:
                parents.extend(attrs["Parent"].split(","))
    missing = [parent for parent in parents if parent not in ids]
    if missing:
        return fail(f"{path}: missing Parent target(s): {', '.join(sorted(set(missing)))}")
    return 0


def validate_fastq_gz(path):
    with gzip.open(path, "rt") as handle:
        lines = [line.rstrip("\n") for line in handle]
    if len(lines) % 4 != 0 or not lines:
        return fail(f"{path}: FASTQ does not contain complete records")
    for offset in range(0, len(lines), 4):
        name, sequence, plus, quality = lines[offset : offset + 4]
        if not name.startswith("@") or plus != "+" or len(sequence) != len(quality):
            return fail(f"{path}: invalid FASTQ record starting at line {offset + 1}")
    return 0


def validate_samplesheets():
    with (VALID / "rnaseq_samplesheet.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    layouts = {row["library_layout"] for row in rows}
    if layouts != {"single", "paired", "long"}:
        return fail("rnaseq_samplesheet.csv must contain single, paired and long layouts")

    with (VALID / "protein_samplesheet.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        protein_path = ROOT / row["filename"]
        if not protein_path.exists():
            return fail(f"protein samplesheet path does not exist: {row['filename']}")
    return 0


def validate_checksums():
    with CHECKSUMS.open() as handle:
        for line in handle:
            expected, relative_path = line.strip().split(maxsplit=1)
            path = ROOT / relative_path
            observed = hashlib.sha256(path.read_bytes()).hexdigest()
            if observed != expected:
                return fail(f"checksum mismatch for {relative_path}")
    return 0


def main():
    for relative_path in REQUIRED_VALID_FILES:
        if not (VALID / relative_path).exists():
            return fail(f"missing fixture: valid/{relative_path}")

    reference = read_fasta(VALID / "reference.fa")
    target = read_fasta(VALID / "target.fa")
    if not reference or not target:
        return fail("reference.fa and target.fa must both contain sequences")

    for gff in ["reference.gff3", "evidence/liftoff_previous_annotations.gff3", "evidence/augustus.hints.gff3"]:
        result = validate_gff3(VALID / gff, target if gff.startswith("evidence/") else reference)
        if result:
            return result

    for fastq in (VALID / "rnaseq").glob("*.fastq.gz"):
        result = validate_fastq_gz(fastq)
        if result:
            return result

    for validator in (validate_samplesheets, validate_checksums):
        result = validator()
        if result:
            return result

    print("Minimal TITAN test data OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
