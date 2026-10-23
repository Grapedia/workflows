#!/usr/bin/env python3
"""Functional enrichment of a gene set against a background (eggNOG COG categories, InterPro entries).

NOT part of the pipeline: an analysis tool. The test set and the background are disjoint gene lists (e.g. genes
new in TITAN vs genes carried over from the previous annotation). For each COG category (a gene with "KL"
counts in K and in L) and each InterPro entry, a two-sided Fisher exact test compares the share of test
genes with the term to the share of background genes; p-values are corrected with Benjamini-Hochberg.
Genes without any eggNOG hit are counted as a separate "no eggNOG annotation" class.

Outputs (in --out-dir): cog_enrichment.tsv, interpro_enrichment.tsv (terms with >= --min-genes test genes).
"""
import argparse
import csv
import math
import os
from collections import Counter, defaultdict

COG = {
    "A": "RNA processing and modification", "B": "Chromatin structure and dynamics",
    "C": "Energy production and conversion", "D": "Cell cycle control, cell division",
    "E": "Amino acid transport and metabolism", "F": "Nucleotide transport and metabolism",
    "G": "Carbohydrate transport and metabolism", "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism", "J": "Translation, ribosomal structure and biogenesis",
    "K": "Transcription", "L": "Replication, recombination and repair",
    "M": "Cell wall/membrane/envelope biogenesis", "N": "Cell motility",
    "O": "Post-translational modification, protein turnover, chaperones",
    "P": "Inorganic ion transport and metabolism", "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "S": "Function unknown", "T": "Signal transduction mechanisms", "U": "Intracellular trafficking and secretion",
    "V": "Defense mechanisms", "W": "Extracellular structures", "Y": "Nuclear structure", "Z": "Cytoskeleton",
}


def log_choose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_two_sided(a, b, c, d):
    """2x2 table [[a, b], [c, d]]: sum of the probabilities <= the observed one."""
    r1, c1, n = a + b, a + c, a + b + c + d
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    def logp(x):
        return log_choose(r1, x) + log_choose(n - r1, c1 - x) - log_choose(n, c1)
    obs = logp(a)
    return min(1.0, sum(math.exp(logp(x)) for x in range(lo, hi + 1) if logp(x) <= obs + 1e-9))


def bh(pvalues):
    order = sorted(range(len(pvalues)), key=lambda i: pvalues[i])
    q, prev = [0.0] * len(pvalues), 1.0
    for rank, i in reversed(list(enumerate(order, 1))):
        prev = min(prev, pvalues[i] * len(pvalues) / rank)
        q[i] = prev
    return q


def read_ids(path):
    return {l.strip() for l in open(path) if l.strip()}


def test_terms(term_genes, test, background, names, min_genes):
    rows = []
    nt, nb = len(test), len(background)
    for term, genes in term_genes.items():
        a, c = len(genes & test), len(genes & background)
        if a < min_genes:
            continue
        p = fisher_two_sided(a, nt - a, c, nb - c)
        ft, fb = a / nt, c / nb
        rows.append([term, names.get(term, ""), a, round(100 * ft, 2), c, round(100 * fb, 2),
                     round(ft / fb, 2) if c else float("inf"), p])
    for row, q in zip(rows, bh([r[-1] for r in rows])):
        row.append(q)
    return sorted(rows, key=lambda r: r[-1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test", required=True, help="gene IDs of the test set")
    ap.add_argument("--background", required=True, help="gene IDs of the background set (disjoint)")
    ap.add_argument("--eggnog", required=True)
    ap.add_argument("--interproscan", required=True)
    ap.add_argument("--min-genes", type=int, default=20)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    test, background = read_ids(args.test), read_ids(args.background)
    universe = test | background

    cog_genes, annotated = defaultdict(set), set()
    for line in open(args.eggnog):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[0] not in universe:
            continue
        annotated.add(f[0])
        for letter in (f[6] if f[6] != "-" else "S"):
            cog_genes[letter].add(f[0])
    cog_genes["none"] = universe - annotated
    names = dict(COG, none="No eggNOG annotation")

    ipr_genes, ipr_names = defaultdict(set), {}
    for f in csv.reader(open(args.interproscan), delimiter="\t"):
        if len(f) > 12 and f[11] not in ("", "-") and f[0] in universe:
            ipr_genes[f[11]].add(f[0])
            ipr_names[f[11]] = f[12]

    os.makedirs(args.out_dir, exist_ok=True)
    header = ["term", "name", "test_genes", "test_pct", "background_genes", "background_pct", "fold", "p", "q"]
    for fname, rows in (("cog_enrichment.tsv", test_terms(cog_genes, test, background, names, 0)),
                        ("interpro_enrichment.tsv", test_terms(ipr_genes, test, background, ipr_names, args.min_genes))):
        with open(os.path.join(args.out_dir, fname), "w") as out:
            out.write("\t".join(header) + "\n")
            for r in rows:
                out.write("\t".join(str(x) for x in r[:7])
                          + f"\t{r[7]:.3g}\t{r[8]:.3g}\n")
    print(f"test {len(test)} genes, background {len(background)} genes")


if __name__ == "__main__":
    main()
