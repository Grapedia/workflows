#!/usr/bin/env python3
import argparse
import json
import re
import sys
from pathlib import Path


DNA_RE = re.compile(r"^[ACGTRYSWKMBDHVNacgtryswkmbdhvn.*-]+$")
PROTEIN_RE = re.compile(r"^[A-Za-z*.-]+$")
VALID_STRANDS = {"+", "-", ".", "?"}
VALID_PHASES = {"0", "1", "2"}


def error(message):
    return f"ERROR: {message}"


def warning(message):
    return f"WARNING: {message}"


def read_fasta(path, sequence_re, label):
    records = {}
    current = None
    chunks = []
    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    records[current] = "".join(chunks)
                current = line[1:].split()[0]
                if not current:
                    raise ValueError(f"{path}: line {line_number}: {label} FASTA header has no sequence id")
                if current in records:
                    raise ValueError(f"{path}: duplicate {label} FASTA sequence id: {current}")
                chunks = []
                continue
            if current is None:
                raise ValueError(f"{path}: line {line_number}: sequence appears before first FASTA header")
            chunks.append(line)
    if current is not None:
        records[current] = "".join(chunks)
    if not records:
        raise ValueError(f"{path}: {label} FASTA contains no sequences")
    empty = [name for name, seq in records.items() if not seq]
    if empty:
        raise ValueError(f"{path}: empty {label} FASTA sequence(s): {', '.join(empty)}")
    invalid = [name for name, seq in records.items() if not sequence_re.match(seq)]
    if invalid:
        raise ValueError(f"{path}: invalid characters in {label} FASTA sequence(s): {', '.join(invalid)}")
    return records


def parse_attributes(raw):
    attrs = {}
    if raw == ".":
        return attrs
    for item in raw.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = value
    return attrs


def validate_gff3(path, genome):
    errors = []
    warnings = []
    ids = {}
    parents = []
    features = []
    feature_counts = {}

    with path.open() as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                errors.append(error(f"{path}: line {line_number}: expected 9 GFF3 columns, found {len(fields)}"))
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            feature_counts[feature_type] = feature_counts.get(feature_type, 0) + 1
            if seqid not in genome:
                errors.append(error(f"{path}: line {line_number}: seqid '{seqid}' is absent from genome FASTA"))
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                errors.append(error(f"{path}: line {line_number}: start/end must be integers"))
                continue
            if start_i < 1 or end_i < start_i or end_i > len(genome[seqid]):
                errors.append(error(f"{path}: line {line_number}: invalid coordinates {start}-{end} for {seqid}"))
            if strand not in VALID_STRANDS:
                errors.append(error(f"{path}: line {line_number}: invalid strand '{strand}'"))
            if feature_type == "CDS" and phase not in VALID_PHASES:
                errors.append(error(f"{path}: line {line_number}: CDS phase must be one of 0, 1 or 2"))
            if feature_type != "CDS" and phase != ".":
                warnings.append(warning(f"{path}: line {line_number}: non-CDS feature has non-dot phase '{phase}'"))
            attrs = parse_attributes(attributes)
            feature_id = attrs.get("ID")
            if feature_id:
                if feature_id in ids:
                    errors.append(error(f"{path}: line {line_number}: duplicate GFF3 ID '{feature_id}'"))
                ids[feature_id] = {
                    "line": line_number,
                    "seqid": seqid,
                    "start": start_i,
                    "end": end_i,
                    "type": feature_type,
                }
            for parent in attrs.get("Parent", "").split(","):
                if parent:
                    parents.append((parent, line_number, seqid, start_i, end_i, feature_type))
            features.append((seqid, start_i, end_i, feature_type, source))

    if not features:
        errors.append(error(f"{path}: GFF3 contains no features"))

    for parent, line_number, seqid, start_i, end_i, feature_type in parents:
        parent_record = ids.get(parent)
        if parent_record is None:
            errors.append(error(f"{path}: line {line_number}: missing Parent target '{parent}'"))
            continue
        if parent_record["seqid"] != seqid:
            errors.append(error(f"{path}: line {line_number}: Parent '{parent}' is on a different seqid"))
            continue
        if start_i < parent_record["start"] or end_i > parent_record["end"]:
            errors.append(error(f"{path}: line {line_number}: {feature_type} coordinates are outside Parent '{parent}'"))

    if feature_counts.get("gene", 0) == 0:
        warnings.append(warning(f"{path}: final annotation contains no gene features"))
    if feature_counts.get("CDS", 0) == 0:
        warnings.append(warning(f"{path}: final annotation contains no CDS features"))

    return errors, warnings, {
        "feature_count": len(features),
        "feature_counts": feature_counts,
        "gff3_ids": len(ids),
    }


def validate_proteins(path, label):
    records = read_fasta(path, PROTEIN_RE, label)
    stop_internal = [
        name for name, seq in records.items()
        if "*" in seq[:-1]
    ]
    warnings = []
    if stop_internal:
        warnings.append(warning(f"{path}: internal stop codon(s) in {label} protein(s): {', '.join(stop_internal)}"))
    return records, warnings


def validate(args):
    errors = []
    warnings = []
    for path, label in [
        (args.genome, "genome FASTA"),
        (args.annotation, "annotation GFF3"),
        (args.proteins_all, "all protein FASTA"),
        (args.proteins_main, "main protein FASTA"),
    ]:
        if not path.exists() or not path.is_file():
            errors.append(error(f"{label} does not exist or is not a file: {path}"))
    if errors:
        return errors, warnings, {}

    try:
        genome = read_fasta(args.genome, DNA_RE, "genome")
    except Exception as exc:
        errors.append(error(str(exc)))
        genome = {}

    if genome:
        gff_errors, gff_warnings, gff_summary = validate_gff3(args.annotation, genome)
        errors.extend(gff_errors)
        warnings.extend(gff_warnings)
    else:
        gff_summary = {}

    try:
        proteins_all, all_warnings = validate_proteins(args.proteins_all, "all")
        warnings.extend(all_warnings)
    except Exception as exc:
        errors.append(error(str(exc)))
        proteins_all = {}

    try:
        proteins_main, main_warnings = validate_proteins(args.proteins_main, "main")
        warnings.extend(main_warnings)
    except Exception as exc:
        errors.append(error(str(exc)))
        proteins_main = {}

    overlap = sorted(set(proteins_all) & set(proteins_main))
    if proteins_all and proteins_main and not overlap:
        warnings.append(warning("main and all protein FASTA files have no shared sequence identifiers"))

    summary = {
        "genome_sequences": len(genome),
        "genome_bases": sum(len(seq) for seq in genome.values()),
        "all_proteins": len(proteins_all),
        "main_proteins": len(proteins_main),
        "shared_protein_ids": len(overlap),
        **gff_summary,
    }
    return errors, warnings, summary


def write_reports(args, errors, warnings, summary):
    status = "fail" if errors else "pass"
    report = {
        "schema_version": "titan.final_annotation_validation.v1",
        "status": status,
        "inputs": {
            "genome": str(args.genome),
            "annotation": str(args.annotation),
            "proteins_all": str(args.proteins_all),
            "proteins_main": str(args.proteins_main),
        },
        "summary": summary,
        "errors": errors,
        "warnings": warnings,
    }
    args.json_report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lines = [
        f"TITAN final annotation validation: {status.upper()}",
        "",
        "Summary:",
    ]
    for key, value in sorted(summary.items()):
        lines.append(f"  {key}: {value}")
    if warnings:
        lines.extend(["", "Warnings:", *[f"  {item}" for item in warnings]])
    if errors:
        lines.extend(["", "Errors:", *[f"  {item}" for item in errors]])
    args.text_report.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Validate TITAN final annotation outputs.")
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--proteins-all", required=True, type=Path)
    parser.add_argument("--proteins-main", required=True, type=Path)
    parser.add_argument("--json-report", default=Path("final_annotation_validation.json"), type=Path)
    parser.add_argument("--text-report", default=Path("final_annotation_validation.txt"), type=Path)
    args = parser.parse_args()

    errors, warnings, summary = validate(args)
    write_reports(args, errors, warnings, summary)
    if errors:
        print("TITAN final annotation validation failed:", file=sys.stderr)
        for item in errors:
            print(f"  {item}", file=sys.stderr)
        return 1
    print("TITAN final annotation validation OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
