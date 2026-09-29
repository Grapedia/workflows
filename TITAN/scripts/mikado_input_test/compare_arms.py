#!/usr/bin/env python3
"""Compare the two arms of the Mikado BRAKER3-input test against the production annotation.

Inputs
  --raw-dir / --tsebra-dir   output dirs of mikado_input_test.nf (data/mikado_input_test/<arm>)
  --prod-mikado              production 04_evidence/mikado/final_mikado_annotation.gff3
  --classes                  data/unsupported_genes_analysis/classes/gene_classes.tsv
  --prod-gff3                production 01_final_annotation/primary/final_annotation.gff3
  --out-dir                  where to write the tables

What is reported
  1. Control: does the `raw` arm reproduce the production Mikado output (gene count, identical loci)?
     If not, the harness differs from production and the test arm cannot be interpreted.
  2. Fate of every production protein-coding gene in each arm (same-strand CDS overlap >= 50 %),
     by evidence class: the test is a success when it loses the `no_signal` /
     `te_overlap_unsupported` genes and (almost) none of the expressed, BUSCO or
     conserved_or_functional genes.
  3. Genes present in an arm but overlapping no production gene ("new").
  4. Completeness and structure: BUSCO summary and AGAT statistics of each arm.
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile

import pandas as pd

CLASS_ORDER = ["expressed", "conserved_or_functional", "te_protein", "te_overlap_unsupported",
               "liftoff_only", "junction_only", "no_signal"]


def gff_genes(path):
    """gene id -> (chrom, start, end, strand); mRNA id -> gene id; CDS intervals per gene."""
    genes, m2g, cds = {}, {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "gene":
                genes[a["ID"]] = (f[0], int(f[3]), int(f[4]), f[6])
            elif f[2] in ("mRNA", "transcript"):
                m2g[a["ID"]] = a["Parent"]
            elif f[2] == "CDS":
                for p in a["Parent"].split(","):
                    g = m2g.get(p)
                    if g:
                        cds.setdefault(g, []).append((f[0], int(f[3]) - 1, int(f[4]), f[6]))
    return genes, cds


def cds_bed(cds, path, strip_id=None):
    with open(path, "w") as out:
        for g, ivs in cds.items():
            for c, s, e, st in ivs:
                out.write("%s\t%d\t%d\t%s\t.\t%s\n" % (c, s, e, g, st))
    subprocess.check_call("sort -k1,1 -k2,2n -o %s %s" % (path, path), shell=True)


def cds_length(cds):
    return {g: sum(e - s for _, s, e, _ in ivs) for g, ivs in cds.items()}


def overlap_fraction(a_cds, b_cds, tmp, tag):
    """For every gene of A: fraction of its CDS bp covered by CDS of genes of B (same strand)."""
    pa, pb = os.path.join(tmp, tag + "_a.bed"), os.path.join(tmp, tag + "_b.bed")
    cds_bed(a_cds, pa)
    cds_bed(b_cds, pb)
    # merge B first so overlapping isoforms are not counted twice
    merged = os.path.join(tmp, tag + "_bm.bed")
    subprocess.check_call(
        "bedtools merge -s -c 6 -o distinct -i %s 2>/dev/null | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t.\\t.\\t\"$4}' > %s"
        % (pb, merged), shell=True)
    out = subprocess.check_output(
        "bedtools intersect -s -a %s -b %s 2>/dev/null | awk -F'\\t' '{d[$4]+=$3-$2} END{for(k in d) print k\"\\t\"d[k]}'"
        % (pa, merged), shell=True).decode()
    ov = {l.split("\t")[0]: int(l.split("\t")[1]) for l in out.strip().split("\n") if l}
    # A's own isoforms overlap each other; cap at 1 after normalising by the union length
    union = {}
    for g, ivs in a_cds.items():
        tot = 0
        last = {}
        for c, s, e, st in sorted(ivs):
            key = (c, st)
            ls, le = last.get(key, (None, None))
            if ls is None or s > le:
                tot += e - s
                last[key] = (s, e)
            elif e > le:
                tot += e - le
                last[key] = (ls, e)
        union[g] = tot
    return pd.Series({g: min(1.0, ov.get(g, 0) / union[g]) if union.get(g) else 0.0 for g in a_cds})


def parse_busco(path):
    if not os.path.exists(path):
        return {}
    txt = open(path).read()
    m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)", txt)
    return dict(zip(["complete", "single", "duplicated", "fragmented", "missing", "n"], m.groups())) if m else {}


def parse_agat(path):
    keys = ["Number of gene", "Number of mrna", "Number of exon", "Number of single exon gene",
            "Number of gene overlapping", "mean mrnas per gene", "mean exons per mrna", "mean gene length (bp)"]
    res = {}
    if not os.path.exists(path):
        return res
    for line in open(path):
        for k in keys:
            # agat_stats.txt has two sections (all isoforms, then longest isoform only): keep the first
            if line.strip().lower().startswith(k.lower()) and k not in res:
                res[k] = line.strip().split()[-1]
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw-dir", required=True)
    ap.add_argument("--tsebra-dir", required=True)
    ap.add_argument("--prod-mikado", required=True)
    ap.add_argument("--prod-gff3", required=True)
    ap.add_argument("--classes", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    report = {}

    gff = lambda d: os.path.join(d, "04_evidence", "mikado", "final_mikado_annotation.gff3")
    arms = {"raw": gff(a.raw_dir), "tsebra": gff(a.tsebra_dir)}
    prod_mk_genes, _ = gff_genes(a.prod_mikado)
    prod_genes, prod_cds = gff_genes(a.prod_gff3)
    classes = pd.read_csv(a.classes, sep="\t", index_col=0)
    prod_cds = {g: v for g, v in prod_cds.items() if g in classes.index}

    with tempfile.TemporaryDirectory() as tmp:
        for arm, path in arms.items():
            if not os.path.exists(path):
                print("missing %s: arm '%s' not finished" % (path, arm), file=sys.stderr)
                continue
            genes, cds = gff_genes(path)
            # 1. control
            same = set(genes.values()) & set(prod_mk_genes.values())
            report[arm] = {"genes": len(genes), "genes_with_cds": len(cds),
                           "loci_identical_to_production_mikado": len(same),
                           "production_mikado_genes": len(prod_mk_genes)}
            # 2. fate of production genes
            frac = overlap_fraction(prod_cds, cds, tmp, arm + "_fate")
            c = classes.join(frac.rename("arm_cds_frac"), how="left")
            c["retained"] = c["arm_cds_frac"].fillna(0) >= 0.5
            tab = c.groupby("evidence_class").agg(
                n=("retained", "size"), retained=("retained", "sum")).reindex(CLASS_ORDER)
            tab["lost"] = tab["n"] - tab["retained"]
            tab["pct_lost"] = (100 * tab["lost"] / tab["n"]).round(1)
            if "busco" in c.columns:
                tab["busco_genes"] = c.groupby("evidence_class")["busco"].sum().reindex(CLASS_ORDER)
                tab["busco_lost"] = c[~c["retained"]].groupby("evidence_class")["busco"].sum().reindex(CLASS_ORDER).fillna(0)
            tab.to_csv(os.path.join(a.out_dir, "fate_by_class_%s.tsv" % arm), sep="\t")
            c[["evidence_class", "arm_cds_frac", "retained"]].to_csv(
                os.path.join(a.out_dir, "fate_by_gene_%s.tsv" % arm), sep="\t")
            report[arm]["fate_by_class"] = tab.fillna(0).astype(int, errors="ignore").to_dict("index")
            # 3. new genes
            back = overlap_fraction(cds, prod_cds, tmp, arm + "_new")
            new = back[back < 0.5]
            pd.Series(new.index).to_csv(os.path.join(a.out_dir, "new_genes_%s.txt" % arm), index=False, header=False)
            report[arm]["genes_not_in_production"] = int(len(new))
            # 4. quality
            report[arm]["busco"] = parse_busco(os.path.join(
                a.raw_dir if arm == "raw" else a.tsebra_dir,
                "01_final_annotation", "quality_report", "busco", "busco_short_summary.txt"))
            report[arm]["agat"] = parse_agat(os.path.join(
                a.raw_dir if arm == "raw" else a.tsebra_dir,
                "01_final_annotation", "quality_report", "agat_stats", "agat_stats.txt"))

    with open(os.path.join(a.out_dir, "comparison.json"), "w") as fh:
        json.dump(report, fh, indent=2, default=str)
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "fate_by_class"} for k, v in report.items()},
                     indent=2, default=str))
    for arm in report:
        print("\n== fate of production genes in arm '%s' ==" % arm)
        print(pd.read_csv(os.path.join(a.out_dir, "fate_by_class_%s.tsv" % arm), sep="\t", index_col=0).to_string())


if __name__ == "__main__":
    main()
