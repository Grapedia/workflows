#!/usr/bin/env python3
"""Evidence metrics (TE, expression, function, non-Vitis homologue) of a reference annotation.

Same definitions as build_gene_evidence_table.py (main isoform = longest CDS, TPM >= 0.5 in >= 1 sample,
"strict" function = Diamond2GO GO / eggNOG description, GO, KEGG / informative InterProScan analysis),
computed over the coding genes (genes with an mRNA) of a GFF3 whose proteins were annotated by
Diamond2GO, eggNOG and InterProScan with the gene ID as query name.

Usage: reference_qc_metrics.py --gff3 ref.gff3 --te-gff3 TEanno.gff3 --tpm gene_tpm_matrix.tsv
           --func-dir functional/v43 --out metrics.json [--nonvitis proteins_main.vs_nonvitis.tsv]
"""
import argparse
import collections
import json
import os
import sys
import tempfile

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_gene_evidence_table import MIN_TPM, NON_INFORMATIVE_IPR, parse_attrs, te_overlap  # noqa: E402


def coding_genes(path):
    genes, mrnas = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = parse_attrs(f[8])
            if f[2] == "gene":
                genes[a["ID"]] = dict(chrom=f[0], start=int(f[3]), end=int(f[4]))
            elif f[2] == "mRNA":
                mrnas[a["ID"]] = dict(gene=a["Parent"], cds=[])
            elif f[2] == "CDS":
                for p in a["Parent"].split(","):
                    if p in mrnas:
                        mrnas[p]["cds"].append((int(f[3]), int(f[4])))
    best = {}
    for m in mrnas.values():
        L = sum(e - s + 1 for s, e in m["cds"])
        if m["gene"] not in best or L > best[m["gene"]][0]:
            best[m["gene"]] = (L, sorted(m["cds"]))
    rows = {}
    for gid, (L, cds) in best.items():
        g = genes[gid]
        rows[gid] = dict(g, cds_intervals=cds, cds_len=L, gene_len=g["end"] - g["start"] + 1)
    return rows


def functional(fdir):
    d = pd.read_csv(os.path.join(fdir, "proteins_main.diamond2go.tsv"), sep="\t", dtype=str)
    d.columns = [c.lstrip("#") for c in d.columns]
    e = pd.read_csv(os.path.join(fdir, "proteins_main.emapper.annotations"), sep="\t", comment="#",
                    header=None, dtype=str, low_memory=False)
    e = e.iloc[:, [0, 7, 9, 11]]
    e.columns = ["g", "desc", "go", "kegg"]
    # InterProScan TSV is ~1.7 GB: read only the columns needed
    ip = pd.read_csv(os.path.join(fdir, "proteins_main.tsv"), sep="\t", header=None, usecols=[0, 3],
                     names=["g", "analysis"], dtype=str, low_memory=False)
    ok = lambda s: s.notna() & ~s.isin(["-", ""])
    fn = dict(
        any=set(d["gene_id"]) | set(e["g"]) | set(ip["g"]),
        strict=set(d.loc[ok(d["GO-term"]), "gene_id"]) | set(e.loc[ok(e["desc"]) | ok(e["go"]) | ok(e["kegg"]), "g"])
        | set(ip.loc[~ip["analysis"].isin(NON_INFORMATIVE_IPR), "g"]))
    return fn


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--te-gff3", required=True)
    ap.add_argument("--tpm", required=True)
    ap.add_argument("--func-dir", required=True)
    ap.add_argument("--nonvitis")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    genes = coding_genes(a.gff3)
    df = pd.DataFrame.from_dict(genes, orient="index")
    n = len(df)
    res = dict(coding_genes=n)
    with tempfile.TemporaryDirectory() as tmp:
        te = te_overlap(genes, df, a.te_gff3, tmp)
    for t in (0.5, 0.8):
        res["te_ge_%d" % (t * 100)] = [int((te["te_frac_cds"] >= t).sum()), 100 * (te["te_frac_cds"] >= t).mean()]
    tpm = pd.read_csv(a.tpm, sep="\t", index_col=0)
    mx = tpm.max(axis=1).reindex(df.index)
    res["tpm_missing_genes"] = int(mx.isna().sum())
    expr = (mx.fillna(0) >= MIN_TPM)
    res["expressed"] = [int(expr.sum()), 100 * expr.mean()]
    fn = functional(a.func_dir)
    any_ = df.index.isin(fn["any"])
    strict = df.index.isin(fn["strict"])
    res["function_any"] = [int(any_.sum()), 100 * any_.mean()]
    res["function_strict"] = [int(strict.sum()), 100 * strict.mean()]
    if a.nonvitis:
        cols = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore".split()
        h = pd.read_csv(a.nonvitis, sep="\t", names=cols)
        h["qcov"] = (h.qend - h.qstart + 1) / h.qlen
        ok = set(h.loc[(h.pident >= 35) & (h.qcov >= 0.5), "qseqid"])
        nv = df.index.isin(ok)
        res["nonvitis"] = [int(nv.sum()), 100 * nv.mean()]
        res["no_evidence_no_liftoff"] = [int((~(any_ | expr.values | nv)).sum()), 100 * (~(any_ | expr.values | nv)).mean()]
    json.dump(res, open(a.out, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
