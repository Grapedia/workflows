#!/usr/bin/env python3
"""Characterise the expressed loci missed by an annotation "without explanation" (unannotated_expression.py).

NOT part of the pipeline: an analysis tool, run on the loci.tsv written by unannotated_expression.py. For each
locus whose reason_missed_by_final starts with "unexplained":
  - exons of its transcripts in the StringTie merges (locus "m0:<gene>" = --stranded-gtf, "m1:" = --unstranded-gtf);
  - overlap with the lncRNA and mRNA exons of another annotation (--other-gff3, e.g. v5.1);
  - Mikado fate: RNA-seq transcripts (star_*, long_*, flair_*) of mikado_prepared.gtf on the locus, and
    whether any transcript of mikado.subloci.gff3 lies on it (none: all excluded by the scoring-file
    requirements); distance to the nearest gene of mikado.loci.gff3 (fragment removal acts within ~2 kb);
  - longest ORF (aa) over its transcripts, both strands for unstranded transcripts.

Usage
  explain_unexplained_loci.py --loci loci.tsv --stranded-gtf merged_stranded.gtf --unstranded-gtf merged_unstranded.gtf \
      --other-gff3 PN40024_5.1_on_T2T_ref.gff3 --mikado-prepared mikado_prepared.gtf \
      --mikado-subloci mikado.subloci.gff3 --mikado-loci mikado.loci.gff3 --genome assembly.fasta \
      --out unexplained_loci.tsv
"""
import argparse
import collections
import os
import re
import subprocess
import tempfile

import pandas as pd

COMP = str.maketrans("ACGTN", "TGCAN")
STOPS = {"TAA", "TAG", "TGA"}
RNASEQ = re.compile(r"^(star_|long_|flair_)")


def bedtools(*args):
    return subprocess.run(["bedtools"] + list(args), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                          universal_newlines=True, check=True).stdout


def write_bed(path, rows):
    with open(path, "w") as out:
        for r in sorted(rows):
            out.write("\t".join(map(str, r)) + "\n")


def load_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip().upper())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def longest_orf(seq):
    best = 0
    for frame in range(3):
        start = None
        for i in range(frame, len(seq) - 2, 3):
            codon = seq[i:i + 3]
            if start is None and codon == "ATG":
                start = i
            elif start is not None and codon in STOPS:
                best = max(best, (i - start) // 3)
                start = None
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for opt in ("loci", "stranded-gtf", "unstranded-gtf", "other-gff3", "mikado-prepared", "mikado-subloci",
                "mikado-loci", "genome", "out"):
        ap.add_argument("--" + opt, required=True)
    a = ap.parse_args()
    tmp = tempfile.mkdtemp()

    loci = pd.read_csv(a.loci, sep="\t")
    u = loci[loci["reason_missed_by_final"].fillna("").str.startswith("unexplained")].copy()
    want = collections.defaultdict(set)
    for locus in u["locus"]:
        want[locus[:2]].add(locus[3:])
    exons, strand = collections.defaultdict(list), {}
    for tag, gtf in (("m0", a.stranded_gtf), ("m1", a.unstranded_gtf)):
        with open(gtf) as fh:
            for line in fh:
                f = line.split("\t")
                if len(f) < 9 or f[2] != "exon":
                    continue
                gene = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
                if gene in want[tag]:
                    key = (tag + ":" + gene, re.search(r'transcript_id "([^"]+)"', f[8]).group(1))
                    exons[key].append((f[0], int(f[3]), int(f[4])))
                    strand[key] = f[6]
    ex_bed = os.path.join(tmp, "exons.bed")
    write_bed(ex_bed, {(c, s - 1, e, k[0]) for k, v in exons.items() for c, s, e in v})

    def loci_hit(bfile):
        return {line.split("\t")[3] for line in bedtools("intersect", "-a", ex_bed, "-b", bfile, "-u").splitlines()}

    for typ in ("lncRNA", "mRNA"):
        tx, rows = set(), []
        with open(a.other_gff3) as fh:
            lines = [line.split("\t") for line in fh if not line.startswith("#")]
        tx = {re.search(r"ID=([^;\n]+)", f[8]).group(1) for f in lines if len(f) > 8 and f[2] == typ}
        rows = [(f[0], int(f[3]) - 1, int(f[4])) for f in lines
                if len(f) > 8 and f[2] == "exon" and re.search(r"Parent=([^;\n]+)", f[8]).group(1) in tx]
        path = os.path.join(tmp, typ + ".bed")
        write_bed(path, rows)
        u["other_" + typ] = u["locus"].isin(loci_hit(path))

    prep = os.path.join(tmp, "prepared_rnaseq.bed")
    rows = []
    with open(a.mikado_prepared) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] == "exon":
                t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
                if RNASEQ.match(t):
                    rows.append((f[0], int(f[3]) - 1, int(f[4])))
    write_bed(prep, rows)
    u["rnaseq_in_mikado_prepare"] = u["locus"].isin(loci_hit(prep))

    sub = os.path.join(tmp, "subloci.bed")
    rows = []
    with open(a.mikado_subloci) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] in ("mRNA", "transcript", "ncRNA"):
                rows.append((f[0], int(f[3]) - 1, int(f[4]), re.search(r"ID=([^;\n]+)", f[8]).group(1)))
    write_bed(sub, rows)
    sub_hits = collections.defaultdict(set)
    for line in bedtools("intersect", "-a", ex_bed, "-b", sub, "-wa", "-wb").splitlines():
        x = line.split("\t")
        sub_hits[x[3]].add(x[7])
    u["in_mikado_subloci"] = u["locus"].isin(sub_hits)
    u["rnaseq_in_mikado_subloci"] = u["locus"].map(lambda l: any(RNASEQ.match(t) for t in sub_hits.get(l, ())))

    genes = os.path.join(tmp, "loci_genes.bed")
    rows = []
    with open(a.mikado_loci) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) > 8 and f[2] in ("gene", "ncRNA_gene"):
                rows.append((f[0], int(f[3]) - 1, int(f[4])))
    write_bed(genes, rows)
    span = os.path.join(tmp, "span.bed")
    write_bed(span, [(r.chrom, r.start - 1, r.end, r.locus) for r in u.itertuples()])
    dist = {}
    for line in bedtools("closest", "-a", span, "-b", genes, "-d", "-t", "first").splitlines():
        x = line.split("\t")
        dist[x[3]] = int(x[-1])
    u["dist_to_mikado_gene"] = u["locus"].map(dist)

    genome = load_fasta(a.genome)
    orf = collections.defaultdict(int)
    for key, iv in exons.items():
        seq = "".join(genome[c][s - 1:e] for c, s, e in sorted(iv, key=lambda z: z[1]))
        rc = seq.translate(COMP)[::-1]
        seqs = [seq] if strand[key] == "+" else [rc] if strand[key] == "-" else [seq, rc]
        orf[key[0]] = max(orf[key[0]], max(longest_orf(s) for s in seqs))
    u["longest_orf_aa"] = u["locus"].map(orf)
    u.to_csv(a.out, sep="\t", index=False)

    for reason, x in u.groupby("reason_missed_by_final"):
        print("%s: %d loci | median samples %d, max TPM %.1f, exonic bp %d" % (
            reason, len(x), x["samples"].median(), x["max_tpm"].median(), x["exonic_bp"].median()))
        print("  other annotation: lncRNA %d, mRNA only %d, neither %d" % (
            x["other_lncRNA"].sum(), (x["other_mRNA"] & ~x["other_lncRNA"]).sum(),
            (~x["other_mRNA"] & ~x["other_lncRNA"]).sum()))
        print("  RNA-seq transcripts in mikado prepare %d; any transcript in subloci %d (RNA-seq %d); "
              "within 2 kb of a Mikado gene %d" % (
                  x["rnaseq_in_mikado_prepare"].sum(), x["in_mikado_subloci"].sum(),
                  x["rnaseq_in_mikado_subloci"].sum(), (x["dist_to_mikado_gene"] <= 2000).sum()))
        print("  longest ORF median %d aa, >= 100 aa %d, >= 200 aa %d" % (
            x["longest_orf_aa"].median(), (x["longest_orf_aa"] >= 100).sum(), (x["longest_orf_aa"] >= 200).sum()))
    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))
    os.rmdir(tmp)


if __name__ == "__main__":
    main()
