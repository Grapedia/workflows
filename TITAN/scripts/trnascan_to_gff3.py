#!/usr/bin/env python3
"""Convert tRNAscan-SE tabular output to minimal GFF3 tRNA features."""

import argparse
import sys
from pathlib import Path
from urllib.parse import quote


HEADER_PREFIXES = ("Sequence", "Name", "--------")


def encode_attribute(value: str) -> str:
    return quote(value, safe="._:-")


def parse_rows(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped.startswith(HEADER_PREFIXES):
                continue

            columns = stripped.split()
            if len(columns) < 9:
                raise ValueError(f"{path}:{line_number}: expected at least 9 columns, found {len(columns)}")

            seqid = columns[0]
            trna_number = columns[1]
            try:
                begin = int(columns[2])
                end = int(columns[3])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid integer coordinate") from error

            trna_type = columns[4]
            anticodon = columns[5]
            score = columns[8]
            strand = "+" if begin <= end else "-"
            start, stop = sorted((begin, end))
            yield seqid, trna_number, start, stop, trna_type, anticodon, score, strand


def convert(input_path: Path, output) -> int:
    print("##gff-version 3", file=output)
    count = 0
    for seqid, trna_number, start, stop, trna_type, anticodon, score, strand in parse_rows(input_path):
        count += 1
        safe_seqid = encode_attribute(seqid)
        feature_id = f"trnascan.{safe_seqid}.tRNA{encode_attribute(trna_number)}"
        attributes = [
            f"ID={feature_id}",
            f"Name={safe_seqid}.tRNA{encode_attribute(trna_number)}",
            f"product=tRNA-{encode_attribute(trna_type)}",
            f"anticodon={encode_attribute(anticodon)}",
        ]
        print(
            "\t".join(
                [
                    seqid,
                    "tRNAscan-SE",
                    "tRNA",
                    str(start),
                    str(stop),
                    score,
                    strand,
                    ".",
                    ";".join(attributes),
                ]
            ),
            file=output,
        )
    return count


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("trnascan_out", type=Path, help="tRNAscan-SE -o table")
    args = parser.parse_args(argv)

    try:
        convert(args.trnascan_out, sys.stdout)
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
