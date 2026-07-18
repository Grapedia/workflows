#!/usr/bin/env python3
"""Summarize gene-level transcriptomic support from Salmon quant.sf files."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def parse_attributes(raw: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in raw.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        attrs[key] = value
    return attrs


def transcript_gene_map(gff_path: Path) -> tuple[set[str], dict[str, str]]:
    genes: set[str] = set()
    tx_to_gene: dict[str, str] = {}
    child_parent: dict[str, str] = {}

    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            feature_type = fields[2]
            attrs = parse_attributes(fields[8])
            feature_id = attrs.get("ID") or attrs.get("transcript_id") or attrs.get("gene_id")
            parent = attrs.get("Parent") or attrs.get("gene_id")

            if feature_type == "gene":
                gene_id = attrs.get("ID") or attrs.get("gene_id")
                if gene_id:
                    genes.add(gene_id)
                continue

            if feature_id and parent:
                child_parent[feature_id] = parent.split(",")[0]

            if feature_type in {"mRNA", "transcript"}:
                transcript_id = attrs.get("ID") or attrs.get("transcript_id")
                gene_id = attrs.get("Parent") or attrs.get("gene_id")
                if transcript_id and gene_id:
                    gene_id = gene_id.split(",")[0]
                    tx_to_gene[transcript_id] = gene_id
                    genes.add(gene_id)

    for transcript_id, gene_id in list(tx_to_gene.items()):
        while gene_id in child_parent:
            gene_id = child_parent[gene_id]
        tx_to_gene[transcript_id] = gene_id
        genes.add(gene_id)

    return genes, tx_to_gene


def sample_name_from_quant(path: Path) -> str:
    if path.name == "quant.sf":
        return path.parent.name.removesuffix("_quant")
    return path.stem


def read_quant_tpm(quant_path: Path) -> dict[str, float]:
    with quant_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "Name" not in reader.fieldnames or "TPM" not in reader.fieldnames:
            raise ValueError(f"{quant_path} is not a Salmon quant.sf file with Name and TPM columns")
        return {row["Name"]: float(row["TPM"]) for row in reader if row.get("Name")}


def summarize(gff_path: Path, quant_paths: list[Path], min_tpm: float) -> dict:
    genes, tx_to_gene = transcript_gene_map(gff_path)
    gene_tpm_by_sample: dict[str, dict[str, float]] = {gene: {} for gene in sorted(genes)}

    for quant_path in quant_paths:
        sample = sample_name_from_quant(quant_path)
        tx_tpm = read_quant_tpm(quant_path)
        sample_gene_tpm: dict[str, float] = defaultdict(float)
        for transcript_id, tpm in tx_tpm.items():
            gene_id = tx_to_gene.get(transcript_id)
            if gene_id:
                sample_gene_tpm[gene_id] += tpm
                genes.add(gene_id)
        for gene_id in sorted(genes):
            gene_tpm_by_sample.setdefault(gene_id, {})[sample] = round(sample_gene_tpm.get(gene_id, 0.0), 6)

    supported = sorted(
        gene_id
        for gene_id, sample_tpms in gene_tpm_by_sample.items()
        if any(tpm >= min_tpm for tpm in sample_tpms.values())
    )
    unsupported = sorted(set(gene_tpm_by_sample) - set(supported))
    total = len(gene_tpm_by_sample)

    return {
        "min_tpm": min_tpm,
        "samples": sorted({sample for sample_tpms in gene_tpm_by_sample.values() for sample in sample_tpms}),
        "total_genes": total,
        "supported_genes": len(supported),
        "unsupported_genes": len(unsupported),
        "supported_gene_percent": round((len(supported) / total * 100.0) if total else 0.0, 4),
        "unsupported_gene_percent": round((len(unsupported) / total * 100.0) if total else 0.0, 4),
        "unsupported_gene_ids": unsupported,
        "gene_tpm": gene_tpm_by_sample,
    }


def write_multiqc(summary: dict, output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("# id: titan_expression_support\n")
        handle.write("# section_name: 'TITAN expression support'\n")
        handle.write("# description: 'Gene-level transcriptomic support from Salmon TPM on the final AEGIS transcriptome.'\n")
        handle.write("# plot_type: 'table'\n")
        handle.write("Metric\tValue\n")
        handle.write(f"Total genes\t{summary['total_genes']}\n")
        handle.write(f"Supported genes\t{summary['supported_genes']}\n")
        handle.write(f"Unsupported genes\t{summary['unsupported_genes']}\n")
        handle.write(f"Supported genes (%)\t{summary['supported_gene_percent']}\n")
        handle.write(f"Minimum TPM\t{summary['min_tpm']}\n")


def write_tpm_matrix(summary: dict, output_path: Path) -> None:
    samples = summary["samples"]
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("gene_id\t" + "\t".join(samples) + "\n")
        for gene_id, sample_tpms in sorted(summary["gene_tpm"].items()):
            handle.write(gene_id + "\t" + "\t".join(str(sample_tpms.get(sample, 0.0)) for sample in samples) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff", required=True, type=Path)
    parser.add_argument("--quant-dirs", nargs="+", required=True, type=Path)
    parser.add_argument("--min-tpm", default=0.5, type=float)
    parser.add_argument("-o", "--output-json", default=Path("expression_support_summary.json"), type=Path)
    parser.add_argument("--multiqc-tsv", default=Path("expression_support_summary_mqc.tsv"), type=Path)
    parser.add_argument("--tpm-matrix", default=Path("gene_tpm_matrix.tsv"), type=Path)
    args = parser.parse_args()

    quant_paths = [path / "quant.sf" if path.is_dir() else path for path in args.quant_dirs]
    missing = [str(path) for path in quant_paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing quant.sf file(s): " + ", ".join(missing))

    summary = summarize(args.gff, quant_paths, args.min_tpm)
    args.output_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_multiqc(summary, args.multiqc_tsv)
    write_tpm_matrix(summary, args.tpm_matrix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
