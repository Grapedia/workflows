#!/usr/bin/env python3
import argparse
import json
from pathlib import Path


def parse_attrs(raw):
    attrs = {}
    for item in raw.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


def gene_intervals(path):
    intervals = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            attrs = parse_attrs(fields[8])
            intervals.append(
                {
                    "seqid": fields[0],
                    "start": start,
                    "end": end,
                    "id": attrs.get("ID", f"{fields[0]}:{start}-{end}"),
                }
            )
    return intervals


def overlaps(left, right):
    return left["seqid"] == right["seqid"] and left["start"] <= right["end"] and right["start"] <= left["end"]


def compare(aegis_gff3, mikado_gff3):
    aegis = gene_intervals(aegis_gff3)
    mikado = gene_intervals(mikado_gff3)
    overlapping = sum(1 for item in mikado if any(overlaps(item, ref) for ref in aegis))
    novel = len(mikado) - overlapping
    return {
        "aegis_gene_count": len(aegis),
        "mikado_gene_count": len(mikado),
        "mikado_genes_overlapping_aegis": overlapping,
        "mikado_genes_not_overlapping_aegis": novel,
    }


def write_reports(summary, json_path, mqc_path):
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with mqc_path.open("w", encoding="utf-8") as handle:
        handle.write("# id: titan_final_annotation_sources\n")
        handle.write("# section_name: 'TITAN final annotation sources'\n")
        handle.write("# description: 'AEGIS and Mikado final GFF3 source comparison.'\n")
        handle.write("# plot_type: 'table'\n")
        handle.write("Metric\tCount\n")
        for key, value in summary.items():
            handle.write(f"{key}\t{value}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare AEGIS and Mikado final annotation GFF3 gene intervals.")
    parser.add_argument("--aegis-gff3", required=True, type=Path)
    parser.add_argument("--mikado-gff3", required=True, type=Path)
    parser.add_argument("--json-report", default=Path("final_annotation_sources.json"), type=Path)
    parser.add_argument("--multiqc-tsv", default=Path("final_annotation_sources_mqc.tsv"), type=Path)
    args = parser.parse_args(argv)
    write_reports(compare(args.aegis_gff3, args.mikado_gff3), args.json_report, args.multiqc_tsv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
