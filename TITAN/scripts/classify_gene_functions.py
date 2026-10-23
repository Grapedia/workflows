#!/usr/bin/env python3
"""Classify genes as TE-related, organellar, other function or unannotated from their functional annotation.

NOT part of the pipeline: an analysis tool (e.g. genes new in TITAN vs genes carried over from the previous
annotation, before function_enrichment.py). For each gene, the text searched is the InterPro entry names, the
signature descriptions (MobiDBLite, Coils and AntiFam left out) and the eggNOG-mapper Description and
Preferred_name of the main protein. Classes are exclusive and tested in this order:

  TE         a transposon keyword (reverse transcriptase, transposase, integrase, gag, copia / gypsy, MuDR,
             hAT, Helitron, En/Spm, RRT15 / ART2, DUF4283, DUF659, PMZ-type, ...)
  organelle  an eggNOG Preferred_name of a plastid / mitochondrial gene (psa, psb, pet, ndh, rpl, rps, rpo,
             ycf, matK, nad, cox, ...) or a description such as photosystem, NAD(P)H-quinone oxidoreductase
             subunit, ribulose bisphosphate carboxylase large chain, maturase K
  other      any other annotation
  none       no annotation

Keyword matching misses TE proteins without a recognisable domain: the TE share is a lower bound.

Outputs in --out-dir: gene_classes.tsv (gene_id, set, class, chrom) and classes.json (counts per set).
"""
import argparse
import collections
import csv
import json
import os
import re

TE = re.compile(
    r"reverse transcriptase|transposase|transposon|retrotranspos|integrase|ribonuclease h domain|"
    r"rnase h domain|gag[- ]|gag-pre|ulp1 protease|hat[, ].*dimeri|mule |mudr|copia|gypsy|ty1|ty3|helitron|"
    r"en/spm|cacta|pol polyprotein|retrovirus|tnt 1-94|harbinger|rrt15|art2|duf4283|duf659|pmz-type", re.I)
ORGANELLE_NAME = re.compile(
    r"^(psa|psb|pet|ndh|atp[abefhi]|rpl\d|rps\d|rpo[abc]|ycf|matk|ccsa|cema|clpp|accd|infa|nad\d|cox\d|cob|"
    r"ccm|mtt|lhb)", re.I)
ORGANELLE_DESC = re.compile(
    r"photosystem|plastoquinone|cytochrome b6|cytochrome b559|NAD\(P\)H-quinone oxidoreductase subunit|"
    r"NADH dehydrogenase subunit|NADH-ubiquinone oxidoreductase chain|ribulose bisphosphate carboxylase large|"
    r"maturase K|DUF825", re.I)
SKIPPED_ANALYSES = ("MobiDBLite", "Coils", "AntiFam")


def read_eggnog(path):
    eggnog = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        eggnog[f[0]] = f
    return eggnog


def read_interproscan(path):
    entries, signatures = collections.defaultdict(set), collections.defaultdict(set)
    for f in csv.reader(open(path), delimiter="\t"):
        if f[3] in SKIPPED_ANALYSES:
            continue
        if len(f) > 12 and f[12] not in ("", "-"):
            entries[f[0]].add(f[12])
        if f[5] not in ("", "-"):
            signatures[f[0]].add(f[5])
    return entries, signatures


def gene_chromosomes(gff3):
    chrom = {}
    for line in open(gff3):
        f = line.split("\t")
        if len(f) == 9 and f[2] == "gene":
            chrom[re.search(r"ID=([^;\n]+)", f[8]).group(1)] = f[0]
    return chrom


def classify(gene, eggnog, entries, signatures):
    e = eggnog.get(gene)
    texts = list(entries[gene]) + list(signatures[gene]) + ([e[7], e[8]] if e else [])
    blob = " | ".join(t for t in texts if t and t != "-")
    if TE.search(blob):
        return "TE"
    if (e and (ORGANELLE_NAME.search(e[8] if e[8] != "-" else "") or ORGANELLE_DESC.search(e[7]))) \
            or any(ORGANELLE_DESC.search(t) for t in entries[gene]):
        return "organelle"
    if blob:
        return "other"
    return "none"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", action="append", required=True, metavar="NAME=FILE",
                    help="gene set to classify (one gene ID per line); repeatable, e.g. new=new_genes.txt")
    ap.add_argument("--eggnog", required=True, help="eggNOG-mapper .emapper.annotations of the main proteins")
    ap.add_argument("--interproscan", required=True, help="InterProScan TSV of the main proteins")
    ap.add_argument("--gff3", required=True, help="annotation GFF3 (gene chromosomes)")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    eggnog = read_eggnog(args.eggnog)
    entries, signatures = read_interproscan(args.interproscan)
    chrom = gene_chromosomes(args.gff3)

    os.makedirs(args.out_dir, exist_ok=True)
    summary = {}
    with open(os.path.join(args.out_dir, "gene_classes.tsv"), "w") as out:
        out.write("gene_id\tset\tclass\tchrom\n")
        for spec in args.set:
            name, path = spec.split("=", 1)
            ids = open(path).read().split()
            counts = collections.Counter()
            for g in ids:
                c = classify(g, eggnog, entries, signatures)
                counts[c] += 1
                out.write(f"{g}\t{name}\t{c}\t{chrom[g]}\n")
            summary[name] = {"n": len(ids), **counts}
            print(name, len(ids), {k: (v, round(100 * v / len(ids), 1)) for k, v in counts.items()})
    json.dump(summary, open(os.path.join(args.out_dir, "classes.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
