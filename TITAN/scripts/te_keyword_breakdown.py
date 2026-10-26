#!/usr/bin/env python3
"""Which transposon keyword puts each gene in the "TE" class of classify_gene_functions.py.

NOT part of the pipeline: an analysis tool. For each gene set of gene_classes.tsv, the genes classed "TE" are
listed with every keyword of classify_gene_functions.TE found in their annotations, and compared with the exact
TE terms removed by the primary_noTE filter (assets/te_function_terms.tsv, InterPro / Pfam / PANTHER
accessions). Also reports their origin, expression and TE (EDTA) overlap from gene_evidence_table.tsv.

Usage
  te_keyword_breakdown.py --classes gene_classes.tsv --eggnog proteins_main.emapper.annotations \
      --interproscan proteins_main.tsv --evidence-table gene_evidence_table.tsv \
      --te-terms assets/te_function_terms.tsv --out te_keyword_genes.tsv
"""
import argparse
import collections
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classify_gene_functions import TE, read_eggnog, read_interproscan  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for opt in ("classes", "eggnog", "interproscan", "evidence-table", "te-terms", "out"):
        ap.add_argument("--" + opt, required=True)
    a = ap.parse_args()

    eggnog = read_eggnog(a.eggnog)
    entries, signatures = read_interproscan(a.interproscan)
    accessions = collections.defaultdict(set)
    with open(a.interproscan) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            accessions[f[0]].update(x for x in (f[4], f[11] if len(f) > 11 else "") if x not in ("", "-"))
    filter_terms = set(pd.read_csv(a.te_terms, sep="\t", comment="#")["term"].dropna())
    classes = pd.read_csv(a.classes, sep="\t", index_col=0)
    evidence = pd.read_csv(a.evidence_table, sep="\t", index_col=0,
                           usecols=["gene_id", "origin", "expressed", "max_tpm", "te_frac_cds"])

    rows = []
    for gene, row in classes[classes["class"] == "TE"].iterrows():
        e = eggnog.get(gene)
        blob = " | ".join(t for t in list(entries[gene]) + list(signatures[gene]) + ([e[7], e[8]] if e else [])
                          if t and t != "-")
        keywords = sorted({m.group(0).lower().strip() for m in TE.finditer(blob)})
        rows.append(dict(gene_id=gene, set=row["set"], chrom=row["chrom"], keywords=";".join(keywords),
                         n_keywords=len(keywords), filter_term=bool(accessions[gene] & filter_terms)))
    df = pd.DataFrame(rows).set_index("gene_id").join(evidence)
    df.to_csv(a.out, sep="\t")

    for name, x in df.groupby("set"):
        kw = collections.Counter(k for v in x["keywords"] for k in v.split(";"))
        single = x.loc[x["n_keywords"] == 1, "keywords"].value_counts()
        expressed = x["expressed"].astype(str).str.lower().eq("true")
        print("== %s: %d TE genes, %d carry a filter term, %d on chr00" % (
            name, len(x), x["filter_term"].sum(), (x["chrom"] == "chr00").sum()))
        print("  genes per keyword:", kw.most_common(15))
        print("  single keyword:", single.head(8).to_dict(), "; several keywords:", int((x["n_keywords"] > 1).sum()))
        print("  origin:", x["origin"].value_counts().to_dict())
        print("  expressed %.1f %%, median max TPM %.2f, CDS >= 50 %% in TE %.1f %%, >= 80 %% %.1f %%" % (
            100 * expressed.mean(), x["max_tpm"].median(),
            100 * (x["te_frac_cds"] >= 0.5).mean(), 100 * (x["te_frac_cds"] >= 0.8).mean()))


if __name__ == "__main__":
    main()
