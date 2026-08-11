#!/usr/bin/env python3
"""Build an EGAPx parameter YAML from TITAN fastp-trimmed read outputs."""

import argparse
from pathlib import Path
from typing import List, Optional, Tuple


def drop_top_level_block(lines: List[str], key: str) -> List[str]:
    output = []  # type: List[str]
    i = 0
    prefix = f"{key}:"

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if line == line.lstrip() and stripped.startswith(prefix):
            i += 1
            while i < len(lines):
                candidate = lines[i]
                candidate_stripped = candidate.strip()
                if (
                    candidate
                    and candidate == candidate.lstrip()
                    and candidate_stripped
                    and not candidate_stripped.startswith("#")
                ):
                    break
                i += 1
            continue

        output.append(line)
        i += 1

    while output and not output[-1].strip():
        output.pop()

    return output


def read_trimmed_reads_manifest(path: Path) -> List[Tuple[str, str, Path, Optional[Path]]]:
    rows = []  # type: List[Tuple[str, str, Path, Optional[Path]]]
    with path.open(encoding="utf-8") as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 4:
                raise ValueError(f"{path}:{line_number}: expected 4 tab-separated fields")
            sample_id, library_layout, read_1, read_2 = fields
            if library_layout not in {"paired", "single"}:
                raise ValueError(
                    f"{path}:{line_number}: library_layout must be 'paired' or 'single', got {library_layout!r}"
                )
            if not sample_id:
                raise ValueError(f"{path}:{line_number}: sampleID is empty")
            if not read_1:
                raise ValueError(f"{path}:{line_number}: read_1 path is empty")
            if library_layout == "paired" and not read_2:
                raise ValueError(f"{path}:{line_number}: paired read row is missing read_2")
            rows.append((sample_id, library_layout, Path(read_1).resolve(), Path(read_2).resolve() if read_2 else None))

    if not rows:
        raise ValueError(f"{path}: no trimmed short-read rows found")
    return rows


def render_short_reads(rows: List[Tuple[str, str, Path, Optional[Path]]]) -> List[str]:
    rendered = ["short_reads:"]
    for sample_id, library_layout, read_1, read_2 in rows:
        rendered.append(f"  - - {sample_id}")
        rendered.append("    - - " + str(read_1))
        if library_layout == "paired":
            rendered.append("      - " + str(read_2))
    return rendered


def build_paramfile(source_yaml: Path, trimmed_reads_tsv: Path, output_yaml: Path) -> None:
    source_lines = source_yaml.read_text(encoding="utf-8").splitlines()
    base_lines = drop_top_level_block(source_lines, "short_reads")
    rows = read_trimmed_reads_manifest(trimmed_reads_tsv)

    rendered = base_lines + [""] + render_short_reads(rows)
    output_yaml.write_text("\n".join(rendered) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Source EGAPx YAML")
    parser.add_argument("--trimmed-reads", type=Path, required=True, help="TSV: sampleID, layout, read1, read2")
    parser.add_argument("--output", type=Path, required=True, help="Derived EGAPx YAML")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_paramfile(args.input, args.trimmed_reads, args.output)


if __name__ == "__main__":
    main()
