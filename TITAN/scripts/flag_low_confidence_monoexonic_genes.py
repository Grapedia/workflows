#!/usr/bin/env python3
"""Classify single-exon genes in the final annotation by supporting evidence.

Ab initio predictors (BRAKER/Helixer) are known to over-predict single-exon
genes in plant genomes. This cross-references every single-exon gene against
evidence TITAN already computes elsewhere in the pipeline (expression,
functional domains, liftoff conservation, BUSCO orthology, TE overlap) and
reports a confidence tier per gene (supported / unsupported+TE-overlap /
unsupported grey-zone) in a TSV+JSON report.

The original final_annotation.gff3 is never touched. In addition, this
writes a "high confidence" GFF3 + main-proteins FASTA variant with every
unsupported (non-"supported"-tier) single-exon gene removed, so its BUSCO/
AGAT stats can be compared side by side against the full annotation's.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import json
import sys
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
            attrs[key] = value
    return attrs


def parse_gff3(gff_path: Path):
    """Single pass over the GFF3: gene records, exon counts per mRNA, and a
    generic child->parent map so any feature ID (CDS, mRNA...) can be walked
    up to its gene."""
    genes: dict[str, dict] = {}
    mrna_to_gene: dict[str, str] = {}
    child_parent: dict[str, str] = {}
    exon_count_by_mrna: dict[str, int] = defaultdict(int)

    with gff_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, _source, feature_type, start, end, _score, strand, _phase, raw_attrs = fields
            attrs = parse_attributes(raw_attrs)
            feature_id = attrs.get("ID")
            parent = attrs.get("Parent")

            if feature_type == "gene" and feature_id:
                genes[feature_id] = {
                    "seqid": seqid,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                }
                continue

            if feature_id and parent:
                child_parent[feature_id] = parent.split(",")[0]

            if feature_type in {"mRNA", "transcript"} and feature_id and parent:
                mrna_to_gene[feature_id] = parent.split(",")[0]
            elif feature_type == "exon" and parent:
                exon_count_by_mrna[parent.split(",")[0]] += 1

    def resolve_to_gene(feature_id: str) -> str | None:
        seen = set()
        current = feature_id
        while current not in genes:
            if current in seen or current not in child_parent:
                return None
            seen.add(current)
            current = child_parent[current]
        return current

    return genes, mrna_to_gene, exon_count_by_mrna, resolve_to_gene


def find_monoexonic_genes(genes, mrna_to_gene, exon_count_by_mrna) -> set[str]:
    mrnas_by_gene: dict[str, list[str]] = defaultdict(list)
    for mrna_id, gene_id in mrna_to_gene.items():
        mrnas_by_gene[gene_id].append(mrna_id)

    monoexonic = set()
    for gene_id in genes:
        mrnas = mrnas_by_gene.get(gene_id)
        if not mrnas:
            continue
        if all(exon_count_by_mrna.get(mrna_id, 0) <= 1 for mrna_id in mrnas):
            monoexonic.add(gene_id)
    return monoexonic


def load_expression(json_summary_path: Path) -> tuple[set[str], float]:
    summary = json.loads(json_summary_path.read_text(encoding="utf-8"))
    return set(summary.get("unsupported_gene_ids", [])), summary.get("min_tpm", 0.0)


def load_liftoff_conserved(correspondence_tsv: Path) -> set[str]:
    conserved = set()
    with correspondence_tsv.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("decision") == "carried_over_from_liftoff":
                conserved.add(row.get("final_gene_id", ""))
    conserved.discard("")
    return conserved


def load_first_column_ids(tsv_path: Path) -> set[str]:
    """Generic reader for InterProScan/eggNOG/Diamond2GO TSVs: the query
    protein ID is always the first column; formats otherwise differ (no
    header at all, '#'-commented header, or a literal 'query' header)."""
    ids = set()
    with tsv_path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            first_col = line.split("\t", 1)[0].strip()
            if not first_col or first_col.lower() == "query":
                continue
            ids.add(first_col)
    return ids


def load_busco_matches(full_table_path: Path) -> set[str]:
    matches = set()
    with full_table_path.open(newline="", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            sequence = fields[2].strip()
            if sequence:
                matches.add(sequence)
    return matches


class TeIndex:
    """Per-chromosome TE intervals, sorted and merged into non-overlapping
    spans so an overlap query is a single bisect (O(log n)) instead of a
    scan over every TE feature per gene (this pipeline has already hit an
    O(n*m) scan on genome-scale data once this session; not repeating it)."""

    def __init__(self, te_gff3_path: Path):
        raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
        with te_gff3_path.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    continue
                seqid, _source, _type, start, end, *_rest = fields
                try:
                    raw[seqid].append((int(start), int(end)))
                except ValueError:
                    continue

        self._spans: dict[str, list[tuple[int, int]]] = {}
        self._starts: dict[str, list[int]] = {}
        for seqid, spans in raw.items():
            spans.sort()
            merged: list[tuple[int, int]] = []
            for start, end in spans:
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            self._spans[seqid] = merged
            self._starts[seqid] = [s for s, _e in merged]

    def overlaps(self, seqid: str, start: int, end: int) -> bool:
        spans = self._spans.get(seqid)
        starts = self._starts.get(seqid)
        if not spans:
            return False
        # Spans are merged and non-overlapping, so at most the one
        # immediately before the insertion point can reach into [start, end].
        idx = bisect.bisect_right(starts, start)
        if idx > 0 and spans[idx - 1][1] >= start:
            return True
        return idx < len(spans) and spans[idx][0] <= end


def build_resolver_for(protein_ids: set[str], resolve_to_gene) -> set[str]:
    """Map a set of protein-derived IDs (which carry a trailing '.prot' from
    the AEGIS FASTA headers) back to gene IDs via the GFF3 parent chain."""
    genes = set()
    for protein_id in protein_ids:
        feature_id = protein_id[:-5] if protein_id.endswith(".prot") else protein_id
        gene_id = resolve_to_gene(feature_id)
        if gene_id:
            genes.add(gene_id)
    return genes


def filter_gff3(gff_path: Path, genes_to_remove: set[str], resolve_to_gene, output_path: Path) -> None:
    """Write every line whose feature does not resolve (via its Parent
    chain) to a removed gene. A line that can't be resolved to any gene at
    all (e.g. a comment, or malformed) is kept - safer to keep an
    unrecognized line than to silently drop real data."""
    with gff_path.open("r", encoding="utf-8") as src, output_path.open("w", encoding="utf-8") as dst:
        for line in src:
            stripped = line.rstrip("\n")
            if not stripped.strip() or stripped.startswith("#"):
                dst.write(line)
                continue
            fields = stripped.split("\t")
            if len(fields) != 9:
                dst.write(line)
                continue
            attrs = parse_attributes(fields[8])
            feature_id = attrs.get("ID")
            gene_id = resolve_to_gene(feature_id) if feature_id else None
            if gene_id in genes_to_remove:
                continue
            dst.write(line)


def filter_proteins_fasta(fasta_path: Path, genes_to_remove: set[str], resolve_to_gene, output_path: Path) -> None:
    """Drop every protein record whose gene is in genes_to_remove. Protein
    headers carry the CDS ID plus a trailing '.prot' (AEGIS convention)."""
    with fasta_path.open("r", encoding="utf-8") as src, output_path.open("w", encoding="utf-8") as dst:
        skip_current = False
        for line in src:
            if line.startswith(">"):
                protein_id = line[1:].split()[0].strip()
                feature_id = protein_id[:-5] if protein_id.endswith(".prot") else protein_id
                gene_id = resolve_to_gene(feature_id)
                skip_current = gene_id in genes_to_remove
            if not skip_current:
                dst.write(line)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff3", required=True, type=Path)
    parser.add_argument("--expression-json", required=True, type=Path)
    parser.add_argument("--te-gff3", required=True, type=Path)
    parser.add_argument("--interproscan-tsv", required=True, type=Path)
    parser.add_argument("--eggnog-tsv", required=True, type=Path)
    parser.add_argument("--diamond2go-tsv", required=True, type=Path)
    parser.add_argument("--liftoff-correspondence-tsv", required=True, type=Path)
    parser.add_argument("--busco-full-table", required=True, type=Path)
    parser.add_argument("--proteins-main-fasta", required=True, type=Path)
    parser.add_argument("--tsv-report", default=Path("monoexonic_gene_confidence.tsv"), type=Path)
    parser.add_argument("--json-summary", default=Path("monoexonic_gene_confidence_summary.json"), type=Path)
    parser.add_argument("--filtered-gff3", default=Path("final_annotation.high_confidence.gff3"), type=Path)
    parser.add_argument(
        "--filtered-proteins-main",
        default=Path("final_annotation_proteins_main.high_confidence.fasta"),
        type=Path,
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    genes, mrna_to_gene, exon_count_by_mrna, resolve_to_gene = parse_gff3(args.gff3)
    monoexonic_genes = find_monoexonic_genes(genes, mrna_to_gene, exon_count_by_mrna)

    unsupported_by_expression, min_tpm = load_expression(args.expression_json)
    liftoff_conserved = load_liftoff_conserved(args.liftoff_correspondence_tsv)

    functional_hit_genes = set()
    functional_hit_genes |= build_resolver_for(load_first_column_ids(args.interproscan_tsv), resolve_to_gene)
    functional_hit_genes |= build_resolver_for(load_first_column_ids(args.eggnog_tsv), resolve_to_gene)
    functional_hit_genes |= build_resolver_for(load_first_column_ids(args.diamond2go_tsv), resolve_to_gene)

    busco_ortholog_genes = build_resolver_for(load_busco_matches(args.busco_full_table), resolve_to_gene)

    te_index = TeIndex(args.te_gff3)

    rows = []
    tier_counts = defaultdict(int)
    busco_monoexonic_contradictions = []

    for gene_id in sorted(monoexonic_genes):
        record = genes[gene_id]
        expressed = gene_id not in unsupported_by_expression
        has_functional_hit = gene_id in functional_hit_genes
        is_liftoff_conserved = gene_id in liftoff_conserved
        is_busco_ortholog = gene_id in busco_ortholog_genes
        te_overlap = te_index.overlaps(record["seqid"], record["start"], record["end"])

        evidence_count = sum([expressed, has_functional_hit, is_liftoff_conserved, is_busco_ortholog])

        if evidence_count > 0:
            tier = "supported"
        elif te_overlap:
            tier = "unsupported_te_overlap"
        else:
            tier = "unsupported_grey_zone"
        tier_counts[tier] += 1

        if is_busco_ortholog and tier != "supported":
            busco_monoexonic_contradictions.append(gene_id)

        rows.append(
            {
                "gene_id": gene_id,
                "seqid": record["seqid"],
                "start": record["start"],
                "end": record["end"],
                "strand": record["strand"],
                "tier": tier,
                "evidence_count": evidence_count,
                "expressed": expressed,
                "has_functional_domain_hit": has_functional_hit,
                "liftoff_conserved": is_liftoff_conserved,
                "busco_ortholog_match": is_busco_ortholog,
                "te_overlap": te_overlap,
            }
        )

    with args.tsv_report.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else [
            "gene_id", "seqid", "start", "end", "strand", "tier", "evidence_count",
            "expressed", "has_functional_domain_hit", "liftoff_conserved",
            "busco_ortholog_match", "te_overlap",
        ], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    genes_to_remove = {row["gene_id"] for row in rows if row["tier"] != "supported"}
    filter_gff3(args.gff3, genes_to_remove, resolve_to_gene, args.filtered_gff3)
    filter_proteins_fasta(args.proteins_main_fasta, genes_to_remove, resolve_to_gene, args.filtered_proteins_main)

    # AGAT's own gene/mRNA stats (what this report is meant to be read
    # alongside) only count protein-coding genes, i.e. genes with at least
    # one mRNA child; `genes` also includes ncRNA-only gene records that
    # have no notion of being mono/multi-exonic in this sense.
    total_genes = len(set(mrna_to_gene.values()))
    total_monoexonic = len(monoexonic_genes)
    summary = {
        "schema_version": "titan.monoexonic_gene_confidence.v1",
        "min_tpm": min_tpm,
        "total_genes": total_genes,
        "monoexonic_genes": total_monoexonic,
        "monoexonic_gene_percent": round((total_monoexonic / total_genes * 100.0) if total_genes else 0.0, 4),
        "tier_counts": dict(tier_counts),
        "regression_check": {
            "description": (
                "Monoexonic genes that are BUSCO single-copy/duplicated/fragmented "
                "orthologs must always land in the 'supported' tier, since a BUSCO "
                "match is itself a positive evidence signal. A non-empty list here "
                "means the evidence mapping has a bug, not that filtering is unsafe."
            ),
            "busco_ortholog_monoexonic_genes": len(
                busco_ortholog_genes & monoexonic_genes
            ),
            "contradictions": busco_monoexonic_contradictions,
        },
        "filtering": {
            "description": (
                "final_annotation.gff3 is never modified. genes_removed lists "
                "every single-exon gene in a non-'supported' tier, dropped (with "
                "all their child features) from --filtered-gff3 / "
                "--filtered-proteins-main so BUSCO/AGAT can be compared against "
                "the full annotation."
            ),
            "genes_removed": len(genes_to_remove),
            "genes_removed_te_overlap": tier_counts.get("unsupported_te_overlap", 0),
            "genes_removed_grey_zone": tier_counts.get("unsupported_grey_zone", 0),
            "filtered_gff3": str(args.filtered_gff3),
            "filtered_proteins_main": str(args.filtered_proteins_main),
        },
        "inputs": {
            "gff3": str(args.gff3),
            "expression_json": str(args.expression_json),
            "te_gff3": str(args.te_gff3),
            "interproscan_tsv": str(args.interproscan_tsv),
            "eggnog_tsv": str(args.eggnog_tsv),
            "diamond2go_tsv": str(args.diamond2go_tsv),
            "liftoff_correspondence_tsv": str(args.liftoff_correspondence_tsv),
            "busco_full_table": str(args.busco_full_table),
            "proteins_main_fasta": str(args.proteins_main_fasta),
        },
    }
    args.json_summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(
        f"{total_monoexonic} single-exon genes out of {total_genes} total "
        f"({summary['monoexonic_gene_percent']}%): "
        f"{tier_counts.get('supported', 0)} supported, "
        f"{tier_counts.get('unsupported_te_overlap', 0)} unsupported+TE-overlap, "
        f"{tier_counts.get('unsupported_grey_zone', 0)} unsupported grey-zone. "
        f"Removed {len(genes_to_remove)} genes into {args.filtered_gff3}.",
        file=sys.stderr,
    )
    if busco_monoexonic_contradictions:
        print(
            f"WARNING: {len(busco_monoexonic_contradictions)} BUSCO-ortholog monoexonic "
            "genes were not classified as supported - check the evidence mapping.",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
