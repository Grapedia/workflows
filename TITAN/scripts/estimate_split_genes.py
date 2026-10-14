#!/usr/bin/env python3
"""Estimate the share of gene models that are split in two (or more) pieces in a TITAN annotation.

NOT part of the pipeline: an analysis tool. It tests every pair of neighbouring genes (same chromosome, same
strand, at most --max-skip genes in between) with three independent kinds of evidence:

  M1 homology   The two genes align to non-overlapping parts of the same reference protein, in protein order
                (the principle of Mender, https://github.com/simrbase/mender). Each gene covers <= 85 % of the
                reference protein, the two alignments meet within +-15 aa, the union covers >= 60 % of it.
                Support = number of distinct reference proteins that show the same tiling.
  M2 long reads Iso-Seq reads whose alignment blocks overlap exons of both genes (primary alignments only).
  M3 expression Pearson correlation of log2(TPM+1) over all RNA-seq samples (both genes expressed in >= 5
                samples), compared with distance-matched control pairs (10-40 genes apart, same strand).

A pair is called "high confidence" when M1 holds (>= 1 reference protein) and M2 or M3 also holds.

Inputs: final GFF3, DIAMOND tables (outfmt 6: qseqid sseqid pident length qlen slen qstart qend sstart send
evalue bitscore, queries = gene IDs), BAM-derived blocks (bedtools bamtobed -split), gene TPM matrix.
Python >= 3.6, pandas, numpy; bedtools on PATH for M2.

Usage
  estimate_split_genes.py --gff3 final_annotation.noTE.gff3 --diamond hits_ref_eu.tsv hits_ref_vs.tsv \
      --long-read-blocks longreads_blocks.bed --tpm-matrix gene_tpm_matrix.tsv --out-dir split_genes
"""
import argparse
import json
import os
import random
import subprocess
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

DIAMOND_COLS = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore".split()


def parse_attrs(field):
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


def read_genes(gff3):
    genes, mrna_gene, exons = {}, {}, defaultdict(list)
    cds_len = defaultdict(int)
    with open(gff3) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            a = parse_attrs(f[8])
            if f[2] == "gene":
                genes[a["ID"]] = dict(chrom=f[0], start=int(f[3]), end=int(f[4]), strand=f[6])
            elif f[2] == "mRNA":
                mrna_gene[a["ID"]] = a["Parent"]
            elif f[2] == "exon":
                for p in a["Parent"].split(","):
                    exons[p].append((int(f[3]), int(f[4])))
            elif f[2] == "CDS":
                for p in a["Parent"].split(","):
                    cds_len[p] += int(f[4]) - int(f[3]) + 1
    main = {}
    for m, g in mrna_gene.items():
        if g not in main or cds_len[m] > cds_len[main[g]]:
            main[g] = m
    return {g: dict(genes[g], exons=sorted(exons[main[g]])) for g in genes if g in main}


def neighbour_pairs(genes, max_skip):
    order = sorted(genes, key=lambda g: (genes[g]["chrom"], genes[g]["start"]))
    pairs = []
    for i, a in enumerate(order):
        for j in range(i + 1, min(i + 2 + max_skip, len(order))):
            b = order[j]
            if genes[b]["chrom"] != genes[a]["chrom"]:
                break
            if genes[b]["strand"] == genes[a]["strand"]:
                pairs.append((a, b, j - i - 1))
    return order, pairs


# ---------------------------------------------------------------- M1
def m1_homology(genes, pairs, diamond_files, gap=15, max_cov_each=0.85, min_union=0.6, min_aln=30):
    hits = defaultdict(dict)   # gene -> {subject: (sstart, send, slen)}
    for path in diamond_files:
        for chunk in pd.read_csv(path, sep="\t", names=DIAMOND_COLS, chunksize=2_000_000):
            chunk = chunk[(chunk["length"] >= min_aln) & chunk["qseqid"].isin(genes)]
            for q, s, ss, se, sl in zip(chunk["qseqid"], chunk["sseqid"], chunk["sstart"], chunk["send"], chunk["slen"]):
                lo, hi = min(ss, se), max(ss, se)
                prev = hits[q].get(s)
                if prev is None or (hi - lo) > (prev[1] - prev[0]):
                    hits[q][s] = (lo, hi, sl)
    support = {}
    for a, b, _skip in pairs:
        ha, hb = hits.get(a), hits.get(b)
        if not ha or not hb:
            continue
        # protein order: + strand = genomic order, - strand = reversed
        first, second = (a, b) if genes[a]["strand"] == "+" else (b, a)
        hf, hs = hits[first], hits[second]
        n = 0
        for s in set(hf) & set(hs):
            (s1, e1, sl), (s2, e2, _) = hf[s], hs[s]
            cov1, cov2 = (e1 - s1 + 1) / sl, (e2 - s2 + 1) / sl
            if cov1 > max_cov_each or cov2 > max_cov_each:
                continue
            if not (-gap <= s2 - e1 - 1 <= gap):
                continue
            union = (max(e1, e2) - min(s1, s2) + 1) / sl
            if union >= min_union:
                n += 1
        if n:
            support[(a, b)] = n
    return support


# ---------------------------------------------------------------- M2
def m2_long_reads(genes, pairs, blocks_bed, out_dir, min_overlap=20):
    exon_bed = os.path.join(out_dir, "gene_exons.bed")
    with open(exon_bed, "w") as out:
        for g, r in genes.items():
            for s, e in r["exons"]:
                out.write("%s\t%d\t%d\t%s\t.\t%s\n" % (r["chrom"], s - 1, e, g, r["strand"]))
    res = subprocess.run(
        "bedtools intersect -wo -a %s -b %s 2>/dev/null" % (blocks_bed, exon_bed),
        shell=True, stdout=subprocess.PIPE, universal_newlines=True, check=True)
    read_genes, same_strand = defaultdict(set), [0, 0]
    ov_by_rg = defaultdict(int)
    for line in res.stdout.splitlines():
        f = line.split("\t")
        read, rstrand, gene, gstrand, ov = f[3], f[5], f[9], f[11], int(f[12])
        ov_by_rg[(read, gene)] += ov
    for (read, gene), ov in ov_by_rg.items():
        if ov >= min_overlap:
            read_genes[read].add(gene)
    n_gene = defaultdict(int)
    for gs in read_genes.values():
        for g in gs:
            n_gene[g] += 1
    span = defaultdict(int)
    for gs in read_genes.values():
        if len(gs) > 1:
            gl = sorted(gs)
            for i in range(len(gl)):
                for j in range(i + 1, len(gl)):
                    span[(gl[i], gl[j])] += 1
    out = {}
    for a, b, _s in pairs:
        k = (a, b) if (a, b) in span else (b, a)
        out[(a, b)] = dict(span=span.get(k, 0), n_a=n_gene.get(a, 0), n_b=n_gene.get(b, 0))
    return out, len(read_genes)


# ---------------------------------------------------------------- M3
def m3_expression(genes, order, pairs, tpm_path, min_samples=5, seed=1, n_control=40000):
    tpm = pd.read_csv(tpm_path, sep="\t", index_col=0)
    tpm = tpm.loc[tpm.index.isin(genes)]
    expr = np.log2(tpm.values + 1.0)
    idx = {g: i for i, g in enumerate(tpm.index)}
    ok = (tpm.values >= 1.0).sum(axis=1) >= min_samples
    z = expr - expr.mean(axis=1, keepdims=True)
    norm = np.sqrt((z ** 2).sum(axis=1))
    norm[norm == 0] = np.nan

    def corr(a, b):
        ia, ib = idx.get(a), idx.get(b)
        if ia is None or ib is None or not (ok[ia] and ok[ib]):
            return np.nan
        return float((z[ia] * z[ib]).sum() / (norm[ia] * norm[ib]))

    adjacent = {(a, b): corr(a, b) for a, b, _s in pairs}
    rng = random.Random(seed)
    control = []
    pos = {g: i for i, g in enumerate(order)}
    for _ in range(n_control):
        a = order[rng.randrange(len(order))]
        j = pos[a] + rng.randint(10, 40)
        if j >= len(order):
            continue
        b = order[j]
        if genes[b]["chrom"] == genes[a]["chrom"] and genes[b]["strand"] == genes[a]["strand"]:
            c = corr(a, b)
            if not np.isnan(c):
                control.append(c)
    return adjacent, np.array(control)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--diamond", nargs="+", required=True)
    ap.add_argument("--long-read-blocks", required=True)
    ap.add_argument("--tpm-matrix", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--max-skip", type=int, default=3, help="genes allowed between the two members of a pair")
    ap.add_argument("--min-span-reads", type=int, default=2)
    ap.add_argument("--min-span-frac", type=float, default=0.10)
    ap.add_argument("--min-corr", type=float, default=0.90)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    genes = read_genes(a.gff3)
    order, pairs = neighbour_pairs(genes, a.max_skip)
    print("genes %d, neighbouring same-strand pairs %d" % (len(genes), len(pairs)), file=sys.stderr)

    m1 = m1_homology(genes, pairs, a.diamond)
    m2, n_reads = m2_long_reads(genes, pairs, a.long_read_blocks, a.out_dir)
    m3, control = m3_expression(genes, order, pairs, a.tpm_matrix)

    rows = []
    for ga, gb, skip in pairs:
        s = m2[(ga, gb)]
        frac = s["span"] / max(min(s["n_a"], s["n_b"]), 1)
        r = m3[(ga, gb)]
        rows.append(dict(
            gene_a=ga, gene_b=gb, genes_between=skip, chrom=genes[ga]["chrom"], strand=genes[ga]["strand"],
            m1_support=m1.get((ga, gb), 0),
            lr_spanning=s["span"], lr_reads_a=s["n_a"], lr_reads_b=s["n_b"], lr_span_frac=round(frac, 3),
            expr_corr="" if np.isnan(r) else round(r, 3),
            m1=m1.get((ga, gb), 0) >= 1,
            m2=s["span"] >= a.min_span_reads and frac >= a.min_span_frac,
            m3=(not np.isnan(r)) and r >= a.min_corr))
    df = pd.DataFrame(rows)
    df["high_confidence"] = df["m1"] & (df["m2"] | df["m3"])
    df.to_csv(os.path.join(a.out_dir, "split_gene_pairs.tsv"), sep="\t", index=False)

    n = len(genes)

    def genes_in(mask):
        sub = df[mask]
        return len(set(sub.gene_a) | set(sub.gene_b))

    adj = df[df.genes_between == 0]
    tested_m3 = adj["expr_corr"].replace("", np.nan).astype(float).notna()
    corr_adj = adj.loc[tested_m3, "expr_corr"].astype(float)
    null_rate = float((control >= a.min_corr).mean()) if len(control) else float("nan")
    summary = dict(
        genes=n, pairs_tested=len(df), pairs_adjacent=len(adj),
        m1_pairs=int(df.m1.sum()), m1_pairs_support_ge3=int((df.m1_support >= 3).sum()),
        m1_genes=genes_in(df.m1), m1_genes_support_ge3=genes_in(df.m1_support >= 3),
        m2_pairs=int(df.m2.sum()), m2_genes=genes_in(df.m2), long_reads_assigned_to_genes=n_reads,
        m2_testable_pairs=int(((df.lr_reads_a >= 3) & (df.lr_reads_b >= 3)).sum()),
        m3_pairs=int(df.m3.sum()), m3_genes=genes_in(df.m3),
        m3_adjacent_tested=int(tested_m3.sum()), m3_adjacent_rate=float((corr_adj >= a.min_corr).mean()),
        m3_control_pairs=int(len(control)), m3_control_rate=null_rate,
        m1_and_m2_pairs=int((df.m1 & df.m2).sum()), m1_and_m3_pairs=int((df.m1 & df.m3).sum()),
        high_conf_pairs=int(df.high_confidence.sum()), high_conf_genes=genes_in(df.high_confidence),
        m1_only_pairs=int((df.m1 & ~df.m2 & ~df.m3).sum()),
        m1_m2_m3_pairs=int((df.m1 & df.m2 & df.m3).sum()),
        params=dict(max_skip=a.max_skip, min_span_reads=a.min_span_reads, min_span_frac=a.min_span_frac, min_corr=a.min_corr),
    )
    with open(os.path.join(a.out_dir, "split_genes_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
