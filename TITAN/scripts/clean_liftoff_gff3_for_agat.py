#!/usr/bin/env python3
"""Clean Liftoff GFF3 records before AGAT sequence extraction."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Gff3Record:
    line: str
    ids: tuple[str, ...]
    parents: tuple[str, ...]


def parse_attributes(raw_attributes: str) -> dict[str, list[str]]:
    attributes: dict[str, list[str]] = {}
    for item in raw_attributes.strip().split(";"):
        if not item:
            continue
        if "=" not in item:
            attributes.setdefault(item, []).append("")
            continue
        key, value = item.split("=", 1)
        attributes.setdefault(key, []).append(value)
    return attributes


def feature_ids(attributes: dict[str, list[str]]) -> tuple[str, ...]:
    return tuple(value for value in attributes.get("ID", []) if value)


def parent_ids(attributes: dict[str, list[str]]) -> tuple[str, ...]:
    parents: list[str] = []
    for value in attributes.get("Parent", []):
        parents.extend(parent for parent in value.split(",") if parent)
    return tuple(parents)


def is_deleted_obsolete_record(line: str, attributes: dict[str, list[str]]) -> bool:
    obsolete = any(value.lower() == "true" for value in attributes.get("obsolete", []))
    return obsolete and "deleted" in line.lower()


def collect_removed_ids(records: list[Gff3Record], seed_removed: set[str]) -> set[str]:
    removed = set(seed_removed)
    changed = True
    while changed:
        changed = False
        for record in records:
            if record.ids and not removed.intersection(record.ids) and removed.intersection(record.parents):
                removed.update(record.ids)
                changed = True
    return removed


def clean_gff3(input_path: Path, output_path: Path, removed_ids_path: Path) -> None:
    comments: list[str] = []
    records: list[Gff3Record] = []
    seed_removed: set[str] = set()

    with input_path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                comments.append(line)
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                continue

            attributes = parse_attributes(fields[8])
            ids = feature_ids(attributes)
            parents = parent_ids(attributes)
            if is_deleted_obsolete_record(line, attributes):
                seed_removed.update(ids)
            records.append(Gff3Record(line=line, ids=ids, parents=parents))

    removed = collect_removed_ids(records, seed_removed)

    with removed_ids_path.open("w", encoding="utf-8") as handle:
        for feature_id in sorted(removed):
            handle.write(f"{feature_id}\n")

    with output_path.open("w", encoding="utf-8") as handle:
        for comment in comments:
            handle.write(f"{comment}\n")
        for record in records:
            if removed.intersection(record.ids) or removed.intersection(record.parents):
                continue
            handle.write(f"{record.line}\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Remove deleted obsolete Liftoff features and descendants from a GFF3 file."
    )
    parser.add_argument("--input", required=True, type=Path, help="Input Liftoff GFF3 file.")
    parser.add_argument("--output", required=True, type=Path, help="Cleaned GFF3 output for AGAT.")
    parser.add_argument(
        "--removed-ids",
        required=True,
        type=Path,
        help="Report of removed feature IDs. The file is created even when no ID is removed.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    clean_gff3(args.input, args.output, args.removed_ids)


if __name__ == "__main__":
    main()
