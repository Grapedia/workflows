#!/usr/bin/env python3
"""Classify overlapping protein-coding gene pairs of a GFF3 by what overlaps.

A pair is counted like AGAT does: the exons of two genes (all isoforms merged) share at
least one base. Each pair is then put in one mutually exclusive class, worst first:

  cds_cds   a CDS base of one gene lies on a CDS base of the other
  cds_utr   a CDS base of one gene lies on an exonic non-CDS (UTR) base of the other
  utr_utr   only UTR bases overlap (no CDS involved)

and split by strand (same / opposite). Genes without an mRNA (ncRNA, lncRNA) are ignored.
"""
import argparse
import json
import sys
from collections import defaultdict


def attrs(raw):
    out = {}
    for item in raw.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1:
            if e > out[-1][1]:
                out[-1][1] = e
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def subtract(a, b):
    """a minus b, both sorted merged interval lists."""
    out = []
    j = 0
    for s, e in a:
        cur = s
        while j < len(b) and b[j][1] < cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] <= e:
            if b[k][0] > cur:
                out.append((cur, b[k][0] - 1))
            cur = max(cur, b[k][1] + 1)
            k += 1
        if cur <= e:
            out.append((cur, e))
    return out


def intersects(a, b):
    i = j = 0
    while i < len(a) and j < len(b):
        if a[i][1] < b[j][0]:
            i += 1
        elif b[j][1] < a[i][0]:
            j += 1
        else:
            return True
    return False


def load(path):
    gene_of = {}
    genes = {}
    mrna_gene = {}
    exons = defaultdict(list)
    cds = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            t = f[2]
            a = attrs(f[8])
            if t == "gene":
                genes[a["ID"]] = (f[0], f[6])
            elif t in ("mRNA", "transcript"):
                if t == "mRNA" and "Parent" in a:
                    mrna_gene[a["ID"]] = a["Parent"].split(",")[0]
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] not in ("exon", "CDS"):
                continue
            a = attrs(f[8])
            for p in a.get("Parent", "").split(","):
                g = mrna_gene.get(p)
                if g is None:
                    continue
                (exons if f[2] == "exon" else cds)[g].append((int(f[3]), int(f[4])))
    coding = {}
    for g in set(mrna_gene.values()):
        if g not in genes:
            continue
        ex = merge(exons[g] + cds[g])
        cd = merge(cds[g])
        if not ex:
            continue
        coding[g] = (genes[g][0], genes[g][1], ex, cd, subtract(ex, cd))
    return coding


def classify(coding):
    by_chr = defaultdict(list)
    for g, (chrom, strand, ex, cd, utr) in coding.items():
        by_chr[chrom].append((ex[0][0], ex[-1][1], g))
    pairs = []
    for chrom, items in by_chr.items():
        items.sort()
        active = []
        for s, e, g in items:
            active = [x for x in active if x[1] >= s]
            for s2, e2, g2 in active:
                if intersects(coding[g][2], coding[g2][2]):
                    pairs.append((g2, g))
            active.append((s, e, g))
    res = []
    for a, b in pairs:
        A, B = coding[a], coding[b]
        if intersects(A[3], B[3]):
            cls = "cds_cds"
        elif intersects(A[3], B[4]) or intersects(B[3], A[4]):
            cls = "cds_utr"
        else:
            cls = "utr_utr"
        res.append((a, b, cls, "same" if A[1] == B[1] else "opposite"))
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--pairs-tsv")
    ap.add_argument("--json")
    args = ap.parse_args()
    coding = load(args.gff)
    res = classify(coding)
    n = len(coding)
    out = {"label": args.label, "coding_genes": n, "pairs": len(res)}
    in_any = set()
    for a, b, c, s in res:
        in_any.update((a, b))
    out["genes_in_overlap"] = len(in_any)
    for cls in ("cds_cds", "cds_utr", "utr_utr"):
        for s in ("same", "opposite"):
            out["%s_%s" % (cls, s)] = sum(1 for r in res if r[2] == cls and r[3] == s)
    # worst class per gene
    rank = {"cds_cds": 0, "cds_utr": 1, "utr_utr": 2}
    worst = {}
    for a, b, c, s in res:
        for g in (a, b):
            if g not in worst or rank[c] < rank[worst[g]]:
                worst[g] = c
    for cls in rank:
        out["genes_worst_" + cls] = sum(1 for v in worst.values() if v == cls)
    cds_same = set()
    for a, b, c, s in res:
        if c == "cds_cds" and s == "same":
            cds_same.update((a, b))
    out["genes_cds_cds_same_strand"] = len(cds_same)
    if args.pairs_tsv:
        with open(args.pairs_tsv, "w") as o:
            o.write("gene_a\tgene_b\tclass\tstrand\n")
            for r in res:
                o.write("\t".join(r) + "\n")
    if args.json:
        with open(args.json, "w") as o:
            json.dump(out, o, indent=1)
    print(json.dumps(out))


if __name__ == "__main__":
    sys.exit(main())
