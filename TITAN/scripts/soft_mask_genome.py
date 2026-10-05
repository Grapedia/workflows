#!/usr/bin/env python3
"""Build a soft-masked genome from the original FASTA and the EDTA hard-masked FASTA.

EDTA (`*.MAKER.masked`) replaces every transposable-element base with `N`. Gene finders such
as BRAKER3/AUGUSTUS/GeneMark-ES expect soft-masking instead (repeats in lower case), which lets
them down-weight repeats without erasing the sequence. This script lower-cases, in the ORIGINAL
sequence, exactly the positions that EDTA turned into `N`.

- Records are matched by FASTA ID (the first word of the header), not by order: EDTA writes the
  records in its own order (chr00 first for PN40024 T2T).
- Both files must contain the same IDs with the same lengths, otherwise the script stops.
- Bases that are `N` in the original assembly stay upper-case `N` (assembly gaps are not repeats).
- Output sequence is wrapped at 80 columns; headers are copied from the original.

Python >= 3.6, standard library only.

Usage: soft_mask_genome.py --genome original.fasta --masked edta.MAKER.masked --out soft.fasta
"""
import argparse
import json
import re
import sys

N_RUN = re.compile(r"[Nn]+")


def read_fasta(path):
    """Yield (header_line_without_>, sequence) one record at a time."""
    header = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif header is not None:
                chunks.append(line.strip())
            elif line.strip():
                raise SystemExit("error: %s: sequence before the first header" % path)
    if header is not None:
        yield header, "".join(chunks)


def record_id(header):
    return header.split()[0] if header.split() else ""


def soft_mask(original, masked):
    """Return (softmasked sequence, number of lower-cased bases)."""
    if len(original) != len(masked):
        raise ValueError("length differs: %d vs %d" % (len(original), len(masked)))
    pieces, last, n_masked = [], 0, 0
    for run in N_RUN.finditer(masked):
        start, end = run.span()
        segment = original[start:end]
        if segment.upper() == "N" * len(segment):
            continue  # native assembly gap, not a repeat
        pieces.append(original[last:start])
        # lower-case the repeat but keep native Ns upper-case
        pieces.append(segment.lower().replace("n", "N"))
        n_masked += end - start - segment.upper().count("N")
        last = end
    pieces.append(original[last:])
    return "".join(pieces), n_masked


def write_record(out, header, seq, width=80):
    out.write(">%s\n" % header)
    for i in range(0, len(seq), width):
        out.write(seq[i:i + width] + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True, help="original (unmasked) FASTA")
    ap.add_argument("--masked", required=True, help="EDTA hard-masked FASTA (repeats = N)")
    ap.add_argument("--out", required=True, help="soft-masked FASTA to write")
    ap.add_argument("--summary", help="optional JSON summary to write")
    args = ap.parse_args()

    masked = {}
    for header, seq in read_fasta(args.masked):
        rid = record_id(header)
        if rid in masked:
            raise SystemExit("error: duplicate ID %s in %s" % (rid, args.masked))
        masked[rid] = seq

    total = lowered = 0
    per_record = {}
    seen = set()
    with open(args.out, "w") as out:
        for header, seq in read_fasta(args.genome):
            rid = record_id(header)
            if rid not in masked:
                raise SystemExit("error: %s is in %s but not in %s" % (rid, args.genome, args.masked))
            try:
                soft, n = soft_mask(seq, masked.pop(rid))
            except ValueError as exc:
                raise SystemExit("error: record %s: %s" % (rid, exc))
            seen.add(rid)
            write_record(out, header, soft)
            total += len(seq)
            lowered += n
            per_record[rid] = round(n / len(seq), 4) if seq else 0.0
    if masked:
        raise SystemExit("error: %d record(s) only in %s: %s"
                         % (len(masked), args.masked, ", ".join(sorted(masked)[:5])))

    frac = lowered / total if total else 0.0
    sys.stderr.write("soft-masked %d of %d bases (%.2f %%) in %d records\n"
                     % (lowered, total, 100 * frac, len(seen)))
    if args.summary:
        with open(args.summary, "w") as fh:
            json.dump({"total_bases": total, "soft_masked_bases": lowered,
                       "soft_masked_fraction": round(frac, 6), "records": len(seen),
                       "fraction_per_record": per_record}, fh, indent=2)


if __name__ == "__main__":
    main()
