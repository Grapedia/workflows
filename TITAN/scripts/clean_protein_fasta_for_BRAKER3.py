#!/usr/bin/env python3
import os
import re
import sys
from pathlib import Path


ALLOWED_AA = set("ABCDEFGHIKLMNPQRSTUVWYZXJO")
GAP_CHARS = set(".-")


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def safe_token(value):
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._-")
    return value or "protein"


def iter_fasta(path):
    header = None
    chunks = []
    header_line = None
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header_line, header, "".join(chunks)
                header = line[1:].strip()
                header_line = line_number
                chunks = []
                if not header:
                    fail(f"{path}: line {line_number}: FASTA header has no sequence id")
            else:
                if header is None:
                    fail(f"{path}: line {line_number}: sequence appears before first FASTA header")
                chunks.append(re.sub(r"\s+", "", line))
        if header is not None:
            yield header_line, header, "".join(chunks)


def clean_sequence(path, header_line, sequence):
    if not sequence:
        fail(f"{path}: record starting line {header_line}: empty protein sequence")

    seq = sequence.upper()
    seq = "".join(aa for aa in seq if aa not in GAP_CHARS)
    seq = seq.rstrip("*")

    if not seq:
        fail(f"{path}: record starting line {header_line}: protein sequence is empty after cleanup")
    if "*" in seq:
        fail(f"{path}: record starting line {header_line}: internal stop codon '*' is not valid for BRAKER3/ProtHint")

    invalid = sorted(set(seq) - ALLOWED_AA)
    if invalid:
        fail(
            f"{path}: record starting line {header_line}: invalid protein character(s): "
            + ", ".join(invalid)
        )
    return seq


if len(sys.argv) != 4:
    print("Usage: clean_protein_fasta_for_BRAKER3.py input.fasta output.fasta prefix_protein_name", file=sys.stderr)
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]
prefix = safe_token(os.path.basename(sys.argv[3]))

input_path = Path(input_file)
records = []
seen_input_ids = set()
duplicate_input_ids = set()
for i, (line_number, header, sequence) in enumerate(iter_fasta(input_path), start=1):
    input_id = header.split()[0]
    if input_id in seen_input_ids:
        duplicate_input_ids.add(input_id)
    seen_input_ids.add(input_id)
    records.append((f"{prefix}_protein{i:06d}", clean_sequence(input_path, line_number, sequence)))

if not records:
    fail(f"{input_path}: FASTA contains no protein sequences")

with open(output_file, "w") as out_f:
    for record_id, sequence in records:
        out_f.write(f">{record_id}\n")
        for start in range(0, len(sequence), 80):
            out_f.write(sequence[start : start + 80] + "\n")

if duplicate_input_ids:
    print(
        "WARNING: duplicate input FASTA IDs were normalized for BRAKER3/ProtHint: "
        + ", ".join(sorted(duplicate_input_ids)[:10]),
        file=sys.stderr,
    )
print(f"Cleaned {len(records)} protein sequence(s) for BRAKER3/ProtHint: {output_file}", file=sys.stderr)
