#!/usr/bin/env python3
"""Compare two Mikado pick outputs (production vs a variant, e.g. pick.clustering.cds_only: true).

NOT part of the pipeline: an analysis tool for scripts/mikado_cds_only_test/run_pick.sh.

Reported for each mikado.loci.gff3
  genes, transcripts; Liftoff (reference) genes with a same-strand CDS overlap;
  recovery of the Liftoff genes lost in production pick (from trace_lost_liftoff_genes.py:
  steps "pick:*", and among them "pick:in_host_utr");
  fusions: primary CDS overlapping >= 2 Liftoff genes; coding genes overlapping another one on the same strand;
  median UTR length of the primary coding transcripts.
Also lists the production genes with no same-strand CDS overlap in the variant (prod_genes_gone.txt).

Usage
  compare_loci.py --prod mikado.loci.gff3 --test variant/mikado.loci.gff3 \
      --liftoff liftoff_previous_annotations.gff3 --lost lost_liftoff_genes.tsv --out-dir compare
"""
import argparse
import collections
import os
import subprocess

import pandas as pd


def attrs(field):
    return dict(kv.split("=", 1) for kv in field.strip().split(";") if "=" in kv)


def read_gff(path):
    genes, m2g, primary = {}, {}, set()
    cds, utr = collections.defaultdict(list), collections.Counter()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = attrs(f[8])
            if f[2] == "gene":
                genes[a["ID"]] = (f[0], int(f[3]) - 1, int(f[4]), f[6])
            elif f[2] in ("mRNA", "transcript", "ncRNA"):
                m2g[a["ID"]] = a["Parent"]
                if a.get("primary", "True") == "True":
                    primary.add(a["ID"])
            elif f[2] == "CDS":
                for p in a["Parent"].split(","):
                    cds[p].append((f[0], int(f[3]) - 1, int(f[4]), f[6]))
            elif "UTR" in f[2]:
                for p in a["Parent"].split(","):
                    utr[p] += int(f[4]) - int(f[3]) + 1
    return genes, m2g, primary, cds, utr


def write_bed(path, rows):
    with open(path, "w") as out:
        for c, s, e, name, st in sorted(rows):
            out.write("%s\t%d\t%d\t%s\t0\t%s\n" % (c, s, e, name, st))


def pairs(bedtools, a, b):
    out = subprocess.run([bedtools, "intersect", "-a", a, "-b", b, "-s", "-wa", "-wb"], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, universal_newlines=True, check=True).stdout
    d = collections.defaultdict(set)
    for line in out.splitlines():
        x = line.split("\t")
        d[x[3]].add(x[9])
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for opt in ("prod", "test", "liftoff", "lost", "out-dir"):
        ap.add_argument("--" + opt, required=True)
    ap.add_argument("--bedtools", default="bedtools")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    _, lo_m2g, _, lo_cds, _ = read_gff(a.liftoff)
    lo_bed = os.path.join(a.out_dir, "liftoff_cds.bed")
    write_bed(lo_bed, [(c, s, e, lo_m2g[t], st) for t, iv in lo_cds.items() for c, s, e, st in iv])
    lo_genes = {lo_m2g[t] for t in lo_cds}
    lost = pd.read_csv(a.lost, sep="\t", index_col=0)
    lost_pick = set(lost.index[lost["step"].str.startswith("pick:")])
    in_utr = set(lost.index[lost["step"] == "pick:in_host_utr"])

    rows, cds_beds = [], {}
    for run, path in (("prod", a.prod), ("test", a.test)):
        genes, m2g, primary, cds, utr = read_gff(path)
        coding = {m2g[t] for t in cds if t in m2g}
        cds_beds[run] = os.path.join(a.out_dir, run + "_cds.bed")
        write_bed(cds_beds[run], [(c, s, e, m2g[t], st) for t, iv in cds.items() if t in m2g for c, s, e, st in iv])
        prim = os.path.join(a.out_dir, run + "_primary_cds.bed")
        write_bed(prim, [(c, s, e, t, st) for t, iv in cds.items() if t in primary for c, s, e, st in iv])
        gb = os.path.join(a.out_dir, run + "_coding_genes.bed")
        write_bed(gb, [(c, s, e, g, st) for g, (c, s, e, st) in genes.items() if g in coding])
        covered = set(pairs(a.bedtools, lo_bed, cds_beds[run])) & lo_genes
        ov = pairs(a.bedtools, gb, gb)
        rows.append(dict(
            run=run, genes=len(genes), coding_genes=len(coding), transcripts=len(m2g),
            liftoff_genes_covered=len(covered), liftoff_genes_total=len(lo_genes),
            recovered_of_lost_in_pick="%d / %d" % (len(covered & lost_pick), len(lost_pick)),
            recovered_of_in_host_utr="%d / %d" % (len(covered & in_utr), len(in_utr)),
            fused_primary_cds_2plus_liftoff=sum(len(v) >= 2 for v in pairs(a.bedtools, prim, lo_bed).values()),
            coding_genes_overlapping_same_strand=sum(len(v - {k}) > 0 for k, v in ov.items()),
            median_primary_utr_nt=float(pd.Series([utr[t] for t in primary if t in cds]).median())))
    df = pd.DataFrame(rows).set_index("run").T
    df.to_csv(os.path.join(a.out_dir, "compare_summary.tsv"), sep="\t")
    print(df.to_string())

    kept = set(pairs(a.bedtools, cds_beds["prod"], cds_beds["test"]))
    _, m2g, _, cds, _ = read_gff(a.prod)
    gone = sorted({m2g[t] for t in cds if t in m2g} - kept)
    with open(os.path.join(a.out_dir, "prod_genes_gone.txt"), "w") as out:
        out.write("\n".join(gone) + "\n")
    print("production genes with no same-strand CDS overlap in the test: %d" % len(gone))


if __name__ == "__main__":
    main()
