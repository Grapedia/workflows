#!/usr/bin/env python3
"""Convert Infernal cmsearch --tblout output to compact Rfam GFF3."""

import argparse
import re
import sys
from pathlib import Path
from urllib.parse import quote


SNRNA_PATTERN = re.compile(r"(^|[^A-Za-z0-9])U([1-6]|11|12)([^A-Za-z0-9]|$)", re.IGNORECASE)


def encode_attribute(value):
    return quote(str(value), safe="._:-")


def infer_rfam_type(query_name, description):
    text = f"{query_name} {description}".lower()
    if "rrna" in text or "ribosomal rna" in text or query_name in {"SSU_rRNA", "LSU_rRNA", "5S_rRNA", "5_8S_rRNA"}:
        return "rRNA"
    if "snrna" in text or SNRNA_PATTERN.search(query_name):
        return "snRNA"
    if "snorna" in text or "snoRNA".lower() in text or "h/aca" in text or "c/d box" in text or "scarna" in text:
        return "snoRNA"
    return "other"


def parse_tblout(path):
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            columns = stripped.split(maxsplit=17)
            if len(columns) < 17:
                raise ValueError(f"{path}:{line_number}: expected at least 17 columns, found {len(columns)}")

            seqid = columns[0]
            query_name = columns[2]
            rfam_id = columns[3]
            try:
                seq_from = int(columns[7])
                seq_to = int(columns[8])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid sequence coordinates") from error
            strand = columns[9]
            score = columns[14]
            evalue = columns[15]
            description = columns[17] if len(columns) > 17 else ""
            start, end = sorted((seq_from, seq_to))
            if strand not in {"+", "-"}:
                strand = "." if seq_from <= seq_to else "-"
            yield seqid, query_name, rfam_id, start, end, strand, score, evalue, description


def convert(input_path, output):
    print("##gff-version 3", file=output)
    count = 0
    seen = {}
    for seqid, query_name, rfam_id, start, end, strand, score, evalue, description in parse_tblout(input_path):
        count += 1
        rfam_type = infer_rfam_type(query_name, description)
        key = (seqid, rfam_id)
        seen[key] = seen.get(key, 0) + 1
        feature_id = f"rfam.{encode_attribute(seqid)}.{encode_attribute(rfam_id)}.{seen[key]}"
        attributes = [
            f"ID={feature_id}",
            f"Name={encode_attribute(query_name)}",
            f"Rfam_ID={encode_attribute(rfam_id)}",
            f"E_value={encode_attribute(evalue)}",
            f"rfam_type={encode_attribute(rfam_type)}",
        ]
        if description:
            attributes.append(f"description={encode_attribute(description)}")
        print(
            "\t".join(
                [
                    seqid,
                    "Infernal/Rfam",
                    rfam_type,
                    str(start),
                    str(end),
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
    parser.add_argument("tblout", type=Path, help="Infernal cmsearch --tblout file")
    args = parser.parse_args(argv)

    try:
        convert(args.tblout, sys.stdout)
    except ValueError as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
