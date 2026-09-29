#!/usr/bin/env python3
"""Add homology / gene-finder-provenance columns to the gene evidence table.

Needs the DIAMOND outputs written by data/unsupported_genes_analysis/diamond/run_diamond.sh
(self.tsv, ref_vs.tsv, ref_eu.tsv) and the TSEBRA-filtered braker.gff3.

Columns added
  tsebra_cds_frac / in_tsebra : fraction of the gene's CDS (same strand) covered by a CDS of
                                braker.gff3 (the TSEBRA-selected BRAKER3 set); in_tsebra >= 0.5.
                                Mikado is fed the RAW AUGUSTUS/GeneMark predictions, not braker.gff3.
  ref_vitis_*   : best hit against Vitales UniProt + Viridiplantae Swiss-Prot (contains Vitis,
                  so it can echo another Vitis annotation)
  ref_nonvitis_*: best hit against non-Vitis eudicot UniProt (independent conservation evidence)
  self_*        : best hit against another gene of this same annotation (paralog / fragment)
  self_hit_supported, self_same_chr, self_dist : properties of that best self-hit
  near_supported_same_strand : a supported same-strand gene lies within 2 kb (split-gene suspect)
"""
import argparse
import os
import subprocess
import sys
import tempfile

import pandas as pd

COLS = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore".split()


def best_hits(path, drop_self=False):
    df = pd.read_csv(path, sep="\t", names=COLS)
    if drop_self:
        df = df[df.qseqid != df.sseqid]
    df = df.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
    df["qcov"] = (df.qend - df.qstart + 1) / df.qlen
    df["scov"] = (df.send - df.sstart + 1) / df.slen
    return df.set_index("qseqid")


def fasta_species(path):
    sp = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                sp[name] = "OS=Vitis " in line
    return sp


def tsebra_overlap(tab, outdir, tmp):
    gff = os.path.join(outdir, "01_final_annotation", "primary", "final_annotation.gff3")
    br = os.path.join(outdir, "04_evidence", "gene_prediction", "braker.gff3")
    fin_bed, br_bed = os.path.join(tmp, "fin.bed"), os.path.join(tmp, "br.bed")
    cmd = (
        r"""awk -F'\t' '$3=="CDS"{split($9,a,"Parent=");sub(/;.*/,"",a[2]);print $1"\t"$4-1"\t"$5"\t"a[2]"\t.\t"$7}' %s"""
        r""" | sort -k1,1 -k2,2n > %s""" % (gff, fin_bed))
    subprocess.check_call(cmd, shell=True)
    cmd = (
        r"""awk -F'\t' '$3=="CDS"{print $1"\t"$4-1"\t"$5"\t.\t.\t"$7}' %s | sort -k1,1 -k2,2n"""
        r""" | bedtools merge -s -c 6 -o distinct -i - 2>/dev/null"""
        r""" | awk '{print $1"\t"$2"\t"$3"\t.\t.\t"$4}' > %s""" % (br, br_bed))
    subprocess.check_call(cmd, shell=True)
    out = subprocess.check_output(
        r"""bedtools intersect -s -a %s -b %s 2>/dev/null | awk -F'\t' '{d[$4]+=$3-$2}END{for(k in d)print k"\t"d[k]}'"""
        % (fin_bed, br_bed), shell=True).decode()
    ov = {}
    for l in out.strip().split("\n"):
        if l:
            k, v = l.split("\t")
            g = k.rsplit("_t", 1)[0]
            ov[g] = max(ov.get(g, 0), int(v))
    frac = pd.Series(ov).reindex(tab.index).fillna(0) / tab["cds_len"].clip(lower=1)
    return frac.clip(upper=1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--diamond-dir", required=True)
    a = ap.parse_args()
    tab = pd.read_csv(a.table, sep="\t", index_col=0)
    drop = [c for c in tab.columns if c.startswith(("ref_", "self_", "tsebra_", "near_")) or c == "in_tsebra"]
    tab = tab.drop(columns=drop)

    with tempfile.TemporaryDirectory() as tmp:
        tab["tsebra_cds_frac"] = tsebra_overlap(tab, a.outdir, tmp)
    tab["in_tsebra"] = tab["tsebra_cds_frac"] >= 0.5

    dd = a.diamond_dir
    vs = best_hits(os.path.join(dd, "ref_vs.tsv"))
    is_vitis = fasta_species(os.path.join(dd, "ref_vitales_swissprot.fasta"))
    vs["is_vitis"] = vs["sseqid"].map(is_vitis)
    tab = tab.join(vs[["pident", "qcov", "bitscore"]].add_prefix("ref_vitis_"))
    tab["ref_vitis_hit_is_vitis"] = tab.index.map(vs["is_vitis"])
    eu_path = os.path.join(dd, "ref_eu.tsv")
    if os.path.exists(eu_path) and os.path.getsize(eu_path) > 0:
        eu = best_hits(eu_path)
        tab = tab.join(eu[["pident", "qcov", "bitscore"]].add_prefix("ref_nonvitis_"))
    else:
        print("WARNING: ref_eu.tsv missing; non-Vitis columns skipped", file=sys.stderr)

    sf = best_hits(os.path.join(dd, "self.tsv"), drop_self=True)
    tab = tab.join(sf[["sseqid", "pident", "qcov", "scov", "bitscore"]].add_prefix("self_"))
    supported = tab["expressed"] | tab["functional_any"] | tab["mapman_any"]
    tab["self_hit_supported"] = tab["self_sseqid"].map(supported)
    tab["self_same_chr"] = tab["self_sseqid"].map(tab["chrom"]) == tab["chrom"]
    tab["self_dist"] = (tab["self_sseqid"].map(tab["start"]) - tab["start"]).abs().where(tab["self_same_chr"])

    # split-gene suspect: supported same-strand neighbour within 2 kb
    t = tab.sort_values(["chrom", "start"])
    g = t.groupby("chrom")
    sup = supported.reindex(t.index)
    prev_sup, next_sup = sup.groupby(t["chrom"]).shift(1), sup.groupby(t["chrom"]).shift(-1)
    prev_str, next_str = g["strand"].shift(1), g["strand"].shift(-1)
    near = (((t["dist_prev"] < 2000) & prev_sup.fillna(False).astype(bool) & (prev_str == t["strand"])) |
            ((t["dist_next"] < 2000) & next_sup.fillna(False).astype(bool) & (next_str == t["strand"])))
    tab["near_supported_same_strand"] = near.reindex(tab.index)
    tab.to_csv(a.table, sep="\t")
    print("updated", a.table, file=sys.stderr)


if __name__ == "__main__":
    main()
