#!/usr/bin/env python3
"""Mark, list and filter genes of a TITAN GFF3 according to their independent evidence.

NOT part of the pipeline: an optional, post-hoc tool. It never modifies the input files and
never renames genes (gene IDs stay stable, so the output can be joined back to any
TITAN table).  Method and validation: docs/user/unsupported_genes.md and the HTML report
docs/reports/annotation_quality_and_gene_count.html.

Input: the per-gene class table written by scripts/classify_unsupported_genes.py
(`gene_classes.tsv`, column `evidence_class`).

Evidence classes and the filter level at which their genes are removed
  expressed, conserved_or_functional   keep      never removed
  no_signal, te_overlap_unsupported    strict    removed by every policy
  liftoff_only, junction_only          broad     removed by `broad` and `broad+te`
  te_protein                           te        removed by `broad+te` only (TE-encoded proteins:
                                                 real proteins, but not host genes)
  (ncRNA genes, or genes absent from the table)  keep

Policies: strict = strict; broad = strict + broad; broad+te = strict + broad + te.

Sub-commands
  mark     write the GFF3 with two attributes added to every gene line:
             titan_evidence_class=<class>;titan_filter_level=<keep|strict|broad|te>
  filter   write a GFF3 without the genes of the chosen policy (and all their mRNA/exon/CDS/UTR
           children); optionally the matching protein FASTA; print a summary
  summary  only print how many genes each policy would remove (no output file)

Examples
  flag_unsupported_genes.py mark   --gff3 final_annotation.gff3 --classes gene_classes.tsv --out marked.gff3
  flag_unsupported_genes.py filter --gff3 final_annotation.gff3 --classes gene_classes.tsv \\
        --policy strict --out final_annotation.strict.gff3 \\
        --proteins final_annotation_proteins_main.fasta --proteins-out proteins_main.strict.fasta
  flag_unsupported_genes.py summary --gff3 final_annotation.gff3 --classes gene_classes.tsv

Python >= 3.6, standard library only.
"""
import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict

LEVEL_OF_CLASS = {
    "expressed": "keep",
    "conserved_or_functional": "keep",
    "no_signal": "strict",
    "te_overlap_unsupported": "strict",
    "liftoff_only": "broad",
    "junction_only": "broad",
    "te_protein": "te",
}
POLICIES = {
    "strict": {"strict"},
    "broad": {"strict", "broad"},
    "broad+te": {"strict", "broad", "te"},
}
CDS_RECORD_SUFFIX = re.compile(r"_t\d+_CDS\d+\.prot$")


def read_classes(path):
    """gene_id -> evidence_class"""
    out = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        first = reader.fieldnames[0]
        if "evidence_class" not in reader.fieldnames:
            raise SystemExit("error: %s has no `evidence_class` column" % path)
        for row in reader:
            out[row[first]] = row["evidence_class"]
    unknown = sorted(set(out.values()) - set(LEVEL_OF_CLASS))
    if unknown:
        raise SystemExit("error: unknown evidence class(es) in %s: %s" % (path, ", ".join(unknown)))
    return out


def level_of(gene_id, classes):
    cls = classes.get(gene_id)
    return cls, LEVEL_OF_CLASS.get(cls, "keep")


def parse_attrs(field):
    attrs = {}
    for part in field.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs


def gff_records(path):
    """Yield (raw_line, fields or None). Comment / blank lines give fields=None."""
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                yield line, None
                continue
            fields = line.rstrip("\n").split("\t")
            yield line, (fields if len(fields) >= 9 else None)


def cmd_mark(a):
    classes = read_classes(a.classes)
    counts = Counter()
    with open(a.out, "w") as out:
        for line, f in gff_records(a.gff3):
            if f is not None and f[2] == "gene":
                gid = parse_attrs(f[8]).get("ID", "")
                cls, level = level_of(gid, classes)
                f[8] = f[8].rstrip(";") + ";titan_evidence_class=%s;titan_filter_level=%s" % (
                    cls or "not_evaluated", level)
                counts[(cls or "not_evaluated", level)] += 1
                out.write("\t".join(f) + "\n")
            else:
                out.write(line)
    print("marked %d genes -> %s" % (sum(counts.values()), a.out), file=sys.stderr)
    for (cls, level), n in sorted(counts.items(), key=lambda x: -x[1]):
        print("  %-26s %-7s %7d" % (cls, level, n), file=sys.stderr)


def genes_to_remove(gff3, classes, policy):
    levels = POLICIES[policy]
    drop, by_class = set(), Counter()
    for _, f in gff_records(gff3):
        if f is not None and f[2] == "gene":
            gid = parse_attrs(f[8]).get("ID", "")
            cls, level = level_of(gid, classes)
            if level in levels:
                drop.add(gid)
                by_class[cls] += 1
    return drop, by_class


def dropped_feature_ids(gff3, drop_genes):
    """IDs of every feature that hangs (at any depth) below a dropped gene."""
    parents = {}
    for _, f in gff_records(gff3):
        if f is None:
            continue
        a = parse_attrs(f[8])
        if "ID" in a and "Parent" in a:
            parents[a["ID"]] = a["Parent"].split(",")
    dropped = set(drop_genes)
    memo = {}

    def below(fid):
        if fid in dropped:
            return True
        if fid in memo:
            return memo[fid]
        memo[fid] = False  # guards against cycles in malformed files
        res = any(below(p) for p in parents.get(fid, ()))
        memo[fid] = res
        return res

    for fid in parents:
        if below(fid):
            dropped.add(fid)
    return dropped


def filter_gff3(a, classes):
    drop_genes, by_class = genes_to_remove(a.gff3, classes, a.policy)
    dropped = dropped_feature_ids(a.gff3, drop_genes)
    kept_genes = removed_lines = 0
    with open(a.out, "w") as out:
        for line, f in gff_records(a.gff3):
            if f is None:
                out.write(line)
                continue
            attrs = parse_attrs(f[8])
            parents = attrs.get("Parent", "").split(",") if "Parent" in attrs else []
            if attrs.get("ID") in dropped or (parents and all(p in dropped for p in parents)):
                removed_lines += 1
                continue
            if f[2] == "gene":
                kept_genes += 1
            out.write(line)
    return drop_genes, by_class, kept_genes, removed_lines


def filter_fasta(fasta_in, fasta_out, drop_genes):
    kept = removed = 0
    keep_record = True
    with open(fasta_in) as fin, open(fasta_out, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                gid = CDS_RECORD_SUFFIX.sub("", line[1:].split()[0])
                keep_record = gid not in drop_genes
                kept += keep_record
                removed += not keep_record
            if keep_record:
                fout.write(line)
    return kept, removed


def cmd_filter(a):
    classes = read_classes(a.classes)
    drop, by_class, kept_genes, removed_lines = filter_gff3(a, classes)
    summary = {
        "policy": a.policy,
        "genes_removed": len(drop),
        "genes_kept": kept_genes,
        "gff3_lines_removed": removed_lines,
        "removed_by_class": dict(by_class),
        "output": a.out,
    }
    if a.proteins:
        if not a.proteins_out:
            raise SystemExit("error: --proteins needs --proteins-out")
        kept, removed = filter_fasta(a.proteins, a.proteins_out, drop)
        summary["proteins_kept"], summary["proteins_removed"] = kept, removed
    if a.removed_list:
        with open(a.removed_list, "w") as fh:
            fh.write("\n".join(sorted(drop)) + ("\n" if drop else ""))
    print(json.dumps(summary, indent=2), file=sys.stderr)


def cmd_summary(a):
    classes = read_classes(a.classes)
    total = sum(1 for _, f in gff_records(a.gff3) if f is not None and f[2] == "gene")
    print("%-10s %8s %8s   %s" % ("policy", "removed", "kept", "removed by class"))
    for policy in POLICIES:
        drop, by_class = genes_to_remove(a.gff3, classes, policy)
        detail = ", ".join("%s=%d" % kv for kv in sorted(by_class.items(), key=lambda x: -x[1]))
        print("%-10s %8d %8d   %s" % (policy, len(drop), total - len(drop), detail))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    sub.required = True  # keyword form needs Python >= 3.7
    for name in ("mark", "filter", "summary"):
        p = sub.add_parser(name)
        p.add_argument("--gff3", required=True, help="final_annotation.gff3")
        p.add_argument("--classes", required=True, help="gene_classes.tsv from classify_unsupported_genes.py")
        if name != "summary":
            p.add_argument("--out", required=True, help="output GFF3")
    f = sub.choices["filter"]
    f.add_argument("--policy", choices=sorted(POLICIES), default="strict")
    f.add_argument("--proteins", help="protein FASTA with CDS-record headers (<gene>_tNNN_CDSN.prot) to filter too")
    f.add_argument("--proteins-out")
    f.add_argument("--removed-list", help="write the removed gene IDs to this file")
    a = ap.parse_args()
    {"mark": cmd_mark, "filter": cmd_filter, "summary": cmd_summary}[a.cmd](a)


if __name__ == "__main__":
    main()
