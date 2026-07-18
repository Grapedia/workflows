#!/usr/bin/env python3
import argparse
from pathlib import Path


def parse_source(raw):
    parts = raw.split(":")
    if len(parts) != 5:
        raise argparse.ArgumentTypeError(
            "--source must be PATH:LABEL:STRANDED:SCORE:IS_REFERENCE"
        )
    path, label, stranded, score, is_reference = parts
    if not path or not label:
        raise argparse.ArgumentTypeError("source PATH and LABEL are required")
    return {
        "path": Path(path),
        "label": label,
        "stranded": stranded,
        "score": score,
        "is_reference": is_reference,
    }


def write_mikado_list(sources, output):
    rows = []
    seen_labels = set()
    for source in sources:
        path = source["path"]
        if not path.exists() or path.stat().st_size == 0:
            continue
        label = source["label"]
        if label in seen_labels:
            raise ValueError(f"duplicate Mikado source label: {label}")
        seen_labels.add(label)
        rows.append(
            [
                str(path),
                label,
                source["stranded"],
                source["score"],
                source["is_reference"],
            ]
        )
    if not rows:
        raise ValueError("no non-empty Mikado input sources")
    output.write_text("\n".join("\t".join(row) for row in rows) + "\n", encoding="utf-8")
    return len(rows)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Create Mikado prepare input list from TITAN evidence files.")
    parser.add_argument("--source", action="append", type=parse_source, required=True)
    parser.add_argument("-o", "--output", type=Path, default=Path("transcript_inputs.tsv"))
    args = parser.parse_args(argv)
    write_mikado_list(args.source, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
