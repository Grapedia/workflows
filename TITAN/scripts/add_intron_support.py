#!/usr/bin/env python3
"""Add exact-intron RNA-seq support to the gene evidence table.

For every multi-exon main isoform of final_annotation.gff3, count how many of its
introns are reproduced *exactly* (same seqid, start, end; strand ignored) by an
intron of an RNA-seq transcript assembly (StringTie / PsiCLASS on STAR alignments,
and StringTie on long reads).  This is a splice-junction check that does not depend
on Salmon quantification, so it can tell "not expressed" apart from "expressed but
multi-mapping / not quantified".

Adds columns: n_introns, n_introns_sr (short-read), n_introns_lr (long-read),
frac_introns_sr, frac_introns_any.

Usage: add_intron_support.py --outdir data/titan_prod_out --table gene_evidence_table.tsv
"""
import argparse
import collections
import glob
import os
import sys

import pandas as pd


def introns_from_gtf(path, into):
    ex = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            tid = f[8].split('transcript_id "', 1)[1].split('"', 1)[0]
            ex[(f[0], tid)].append((int(f[3]), int(f[4])))
    for (chrom, _), es in ex.items():
        if len(es) < 2:
            continue
        es.sort()
        for (s1, e1), (s2, e2) in zip(es, es[1:]):
            into.add((chrom, e1 + 1, s2 - 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--table", required=True, help="gene_evidence_table.tsv (rewritten in place)")
    a = ap.parse_args()
    ta = os.path.join(a.outdir, "04_evidence", "transcript_assemblies")
    sr, lr = set(), set()
    for g in sorted(glob.glob(os.path.join(ta, "merged_star_*.gtf"))):
        print("short-read introns:", os.path.basename(g), file=sys.stderr)
        introns_from_gtf(g, sr)
    for g in sorted(glob.glob(os.path.join(ta, "merged_minimap2_*.gtf"))):
        print("long-read introns:", os.path.basename(g), file=sys.stderr)
        introns_from_gtf(g, lr)
    print("distinct introns: short-read %d, long-read %d" % (len(sr), len(lr)), file=sys.stderr)

    tab = pd.read_csv(a.table, sep="\t", index_col=0)
    main_mrna = tab["main_mrna"].to_dict()
    wanted = set(main_mrna.values())
    ex = collections.defaultdict(list)
    chrom_of = {}
    with open(os.path.join(a.outdir, "01_final_annotation", "primary", "final_annotation.gff3")) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            par = f[8].split("Parent=", 1)[1].split(";", 1)[0]
            if par in wanted:
                ex[par].append((int(f[3]), int(f[4])))
                chrom_of[par] = f[0]
    res = {}
    for gid, mid in main_mrna.items():
        es = sorted(ex.get(mid, []))
        introns = [(chrom_of[mid], e1 + 1, s2 - 1) for (s1, e1), (s2, e2) in zip(es, es[1:])] if es else []
        n = len(introns)
        n_sr = sum(i in sr for i in introns)
        n_lr = sum(i in lr for i in introns)
        n_any = sum((i in sr) or (i in lr) for i in introns)
        res[gid] = dict(n_introns=n, n_introns_sr=n_sr, n_introns_lr=n_lr,
                        frac_introns_sr=(n_sr / n) if n else float("nan"),
                        frac_introns_any=(n_any / n) if n else float("nan"))
    new = pd.DataFrame.from_dict(res, orient="index")
    tab = tab.drop(columns=[c for c in new.columns if c in tab.columns]).join(new)
    tab.to_csv(a.table, sep="\t")
    print("updated", a.table, file=sys.stderr)


if __name__ == "__main__":
    main()
