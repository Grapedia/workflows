#!/usr/bin/env python3
"""Audit the removal of every gene carrying one of the functional terms marked for deletion in a screening sheet.

NOT part of the pipeline: an analysis tool, run on the sheet written by list_frequent_functions.py once the
"à supprimer (oui/non)" column has been filled in (sheet "all"). A gene is removed when its main protein carries
at least one term marked "oui", whatever its other annotations:

  InterPro / Pfam / PANTHER   matched on the accession (InterProScan TSV)
  eggNOG                      matched on the exact Description
  Diamond2GO                  matched on the exact hit description (accession removed)

Reported: removed genes (new vs carried over from the previous annotation, expression, other functions kept),
BUSCO before / after (a BUSCO whose complete copies are all removed becomes missing, a duplicated one left with
a single copy becomes single-copy), genes with a function before / after, and the TE-keyword annotations left on
the kept genes. OMArk is not recomputed (rerun it on the filtered proteins).

Outputs: <out>.removed_genes.tsv (one row per removed gene, with the matched terms) and <out>.summary.txt.
"""
import argparse
import collections
import csv
import os
import sys

import openpyxl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from list_frequent_functions import TE  # noqa: E402

NONINF = ("MobiDBLite", "Coils", "AntiFam")
DECISION_COL = "à supprimer (oui/non)"


def read_selection(xlsx, column):
    ws = openpyxl.load_workbook(xlsx, data_only=True)["all"]
    rows = list(ws.iter_rows(values_only=True))
    col = {h: i for i, h in enumerate(rows[0])}
    selected = collections.defaultdict(dict)
    for r in rows[1:]:
        if str(r[col[column]]).strip().lower() == "oui":
            key = r[col["term"]] if r[col["term"]] else r[col["name"]]
            selected[r[col["source"]]][key] = r[col["name"]]
    return selected


def d2g_description(field):
    s = field.strip()
    return s.split(" ", 1)[-1].strip() if " " in s else s


def gene_ids(gff3):
    out = []
    for line in open(gff3):
        f = line.split("\t")
        if len(f) == 9 and f[2] == "gene":
            out.append(f[8].split("ID=")[1].split(";")[0].strip())
    return out


def busco_table(path):
    """BUSCO id -> (status, set of genes)."""
    table = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 3:
            continue
        status, genes = table.get(f[0], (f[1], set()))
        genes.add(f[2].split("_t")[0])
        table[f[0]] = (status, genes)
    return table


def busco_counts(table, n_total, removed):
    c = collections.Counter()
    for status, genes in table.values():
        left = genes - removed
        if not left:
            continue
        if status == "Fragmented":
            c["F"] += 1
        elif len(left) > 1:
            c["D"] += 1
        else:
            c["S"] += 1
    c["C"] = c["S"] + c["D"]
    c["M"] = n_total - c["C"] - c["F"]
    return {k: round(100 * c[k] / n_total, 1) for k in "CSDFM"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selection", required=True, help="screening .xlsx with the 'à supprimer (oui/non)' column")
    ap.add_argument("--column", default=DECISION_COL, help="decision column of the 'all' sheet (cells 'oui')")
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--eggnog", required=True)
    ap.add_argument("--interproscan", required=True)
    ap.add_argument("--diamond2go", required=True)
    ap.add_argument("--busco-full-table", required=True)
    ap.add_argument("--busco-n", type=int, default=1990, help="number of BUSCOs in the lineage")
    ap.add_argument("--tpm-matrix", required=True, help="gene_tpm_matrix.tsv (expressed = max TPM >= --min-tpm)")
    ap.add_argument("--min-tpm", type=float, default=0.5)
    ap.add_argument("--new-genes", help="IDs of genes absent from the previous annotation")
    ap.add_argument("--out", required=True, help="output prefix")
    args = ap.parse_args()

    sel = read_selection(args.selection, args.column)
    genes = gene_ids(args.gff3)
    ipr, egg, d2g = collections.defaultdict(list), {}, collections.defaultdict(list)
    for f in csv.reader(open(args.interproscan), delimiter="\t"):
        ipr[f[0]].append(f)
    for line in open(args.eggnog):
        if not line.startswith("#"):
            f = line.rstrip("\n").split("\t")
            egg[f[0]] = f
    with open(args.diamond2go) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 1:
                d2g[f[0]].append(d2g_description(f[1]))
    max_tpm = {}
    with open(args.tpm_matrix) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            max_tpm[f[0]] = max(float(x) for x in f[1:]) if len(f) > 1 else 0.0
    new = set(open(args.new_genes).read().split()) if args.new_genes else None
    busco = busco_table(args.busco_full_table)
    busco_genes = {g for _, gs in busco.values() for g in gs}

    def matched(g):
        terms = set()
        for f in ipr[g]:
            if len(f) > 12 and f[11] in sel["InterPro"]:
                terms.add(f"InterPro:{f[11]}")
            if f[3] in ("Pfam", "PANTHER") and f[4] in sel[f[3]]:
                terms.add(f"{f[3]}:{f[4]}")
        if g in egg and egg[g][7] in sel["eggNOG"]:
            terms.add(f"eggNOG:{egg[g][7]}")
        terms.update(f"Diamond2GO:{d}" for d in d2g[g] if d in sel["Diamond2GO"])
        return terms

    def has_function(g):
        return bool([f for f in ipr[g] if f[3] not in NONINF] or g in egg or d2g[g])

    def other_function(g):
        """an annotation that is neither selected nor TE-like."""
        texts = [f[12] if len(f) > 12 and f[11] != "-" else f[5] for f in ipr[g]
                 if f[3] not in NONINF and not (len(f) > 12 and f[11] in sel["InterPro"])
                 and not (f[3] in ("Pfam", "PANTHER") and f[4] in sel[f[3]])]
        if g in egg and egg[g][7] not in sel["eggNOG"]:
            texts.append(egg[g][7])
        texts += [d for d in d2g[g] if d not in sel["Diamond2GO"]]
        return [t for t in texts if t and t != "-" and not TE.search(t)]

    removed = {g: matched(g) for g in genes}
    removed = {g: t for g, t in removed.items() if t}
    rset = set(removed)
    kept = [g for g in genes if g not in rset]

    with open(args.out + ".removed_genes.tsv", "w") as out:
        out.write("gene_id\tnew_gene\tbusco\tmax_tpm\tother_non_te_functions\tmatched_terms\n")
        for g in sorted(removed):
            other = other_function(g)
            out.write("\t".join([g, "" if new is None else str(g in new), str(g in busco_genes),
                                 f"{max_tpm.get(g, 0):.2f}", str(len(other)),
                                 " | ".join(sorted(removed[g]))]) + "\n")

    lines = []
    n = len(genes)
    lines.append(f"Selected terms: " + ", ".join(f"{k} {len(v)}" for k, v in sel.items()))
    lines.append(f"Genes: {n} -> {len(kept)} ({len(rset)} removed, {100 * len(rset) / n:.1f} %)")
    if new is not None:
        lines.append(f"  removed: new {len(rset & new)}, carried over {len(rset - new)}")
    expr = [g for g in rset if max_tpm.get(g, 0) >= args.min_tpm]
    lines.append(f"  removed and expressed (max TPM >= {args.min_tpm}): {len(expr)}")
    lines.append(f"  removed with another non-TE function: {sum(1 for g in rset if other_function(g))}")
    rb = sorted(rset & busco_genes)
    lines.append(f"  removed BUSCO genes: {len(rb)} {rb}")
    b0, b1 = busco_counts(busco, args.busco_n, set()), busco_counts(busco, args.busco_n, rset)
    lines.append("BUSCO before: C:{C}% [S:{S}%, D:{D}%], F:{F}%, M:{M}%".format(**b0))
    lines.append("BUSCO after : C:{C}% [S:{S}%, D:{D}%], F:{F}%, M:{M}%".format(**b1))
    lost = [b for b, (s, gs) in busco.items() if s != "Fragmented" and not gs - rset]
    lines.append(f"  BUSCOs lost: {len(lost)} {sorted(lost)}")
    f0 = sum(has_function(g) for g in genes)
    f1 = sum(has_function(g) for g in kept)
    lines.append(f"Genes with a function: {100 * f0 / n:.1f} % -> {100 * f1 / len(kept):.1f} %")
    e0 = sum(max_tpm.get(g, 0) >= args.min_tpm for g in genes)
    e1 = sum(max_tpm.get(g, 0) >= args.min_tpm for g in kept)
    lines.append(f"Expressed genes: {100 * e0 / n:.1f} % -> {100 * e1 / len(kept):.1f} %")
    left = collections.Counter()
    n_left = 0
    for g in kept:
        texts = [f[12] for f in ipr[g] if len(f) > 12] + [f[5] for f in ipr[g]]
        texts += ([egg[g][7]] if g in egg else []) + d2g[g]
        hits = {TE.search(t).group(0).lower() for t in texts if t and t != "-" and TE.search(t)}
        if hits:
            n_left += 1
            left.update(hits)
    lines.append(f"Kept genes still carrying a TE-keyword annotation: {n_left}; top: {left.most_common(10)}")
    open(args.out + ".summary.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
