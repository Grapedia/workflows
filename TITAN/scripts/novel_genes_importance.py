#!/usr/bin/env python3
"""Show that the genes added by TITAN on top of the reference (Liftoff) annotation are real and useful.

Groups (from gene_evidence_table.tsv, column liftoff_carried):
  reference = genes carried over from the Liftoff annotation (validated reference)
  new       = all other genes of final_annotation.gff3
  new_noTE  = new genes with < 50 % of their CDS in an EDTA TE

For each group: expression summary (max TPM over samples, mean TPM, number of samples expressed),
share of genes expressed / with a functional hit / with a non-TE function, breakdown by origin,
and the functional profile (eggNOG COG categories, top eggNOG descriptions, top KEGG pathways),
as a % of the genes of the group.

Usage:
    novel_genes_importance.py --table data/unsupported_genes_analysis/gene_evidence_table.tsv \
        --outdir data/titan_prod_out --out-dir data/unsupported_genes_analysis/novel_genes
Python >= 3.6, pandas required.  Writes summary.md and the TSV tables next to it.
"""
import argparse
import os
import re

import pandas as pd

MIN_TPM = 0.5
CDS_SUFFIX = re.compile(r"_t\d+_CDS\d+\.prot$")
COG = {
    "J": "Translation, ribosomal structure", "A": "RNA processing and modification",
    "K": "Transcription", "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics", "D": "Cell cycle, division, chromosome partitioning",
    "Y": "Nuclear structure", "V": "Defense mechanisms", "T": "Signal transduction",
    "M": "Cell wall/membrane/envelope biogenesis", "N": "Cell motility", "Z": "Cytoskeleton",
    "W": "Extracellular structures", "U": "Intracellular trafficking and secretion",
    "O": "Post-translational modification, chaperones", "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism", "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism", "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism", "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis", "R": "General function prediction only",
    "S": "Function unknown",
}


def describe(s, label):
    q = s.quantile([0.25, 0.5, 0.75, 0.9])
    return {"group": label, "n": len(s), "min": s.min(), "q1": q[0.25], "median": q[0.5],
            "mean": s.mean(), "q3": q[0.75], "p90": q[0.9], "max": s.max()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True, help="gene_evidence_table.tsv")
    ap.add_argument("--outdir", required=True, help="TITAN output dir (for the eggNOG annotations)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--top", type=int, default=25)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    t = pd.read_csv(a.table, sep="\t", index_col="gene_id")
    t["group"] = t["liftoff_carried"].map({True: "reference", False: "new"})
    groups = {"reference": t[t.liftoff_carried], "new": t[~t.liftoff_carried],
              "new_noTE": t[~t.liftoff_carried & (t.te_frac_cds.fillna(0) < 0.5)]}

    # --- expression
    rows, stats = [], []
    for g, d in groups.items():
        for col in ("max_tpm", "mean_tpm", "n_samples_expr"):
            stats.append(dict(describe(d[col], g), metric=col))
        rows.append({"group": g, "n_genes": len(d),
                     "pct_expressed_tpm%.1f" % MIN_TPM: 100 * d.expressed.mean(),
                     "pct_expressed_>=3_samples": 100 * (d.n_samples_expr >= 3).mean(),
                     "pct_functional_any": 100 * d.functional_any.mean(),
                     "pct_functional_strict": 100 * d.functional_strict.mean(),
                     "pct_expressed_or_functional": 100 * (d.expressed | d.functional_any).mean(),
                     "pct_nonvitis_homolog": 100 * ((d.ref_nonvitis_qcov.fillna(0) >= 0.5) & (d.ref_nonvitis_pident.fillna(0) >= 35)).mean(),
                     "pct_omamer_hit": 100 * d.omamer_hit.astype(bool).mean(),
                     "pct_TE_cds>=50": 100 * (d.te_frac_cds.fillna(0) >= 0.5).mean(),
                     "median_prot_len": d.prot_len.median()})
    pd.DataFrame(stats)[["metric", "group", "n", "min", "q1", "median", "mean", "q3", "p90", "max"]] \
        .to_csv(os.path.join(a.out_dir, "expression_summary.tsv"), sep="\t", index=False, float_format="%.3f")
    support = pd.DataFrame(rows)
    support.to_csv(os.path.join(a.out_dir, "support_summary.tsv"), sep="\t", index=False, float_format="%.2f")
    by_origin = pd.crosstab(t.origin, t.group)
    by_origin.to_csv(os.path.join(a.out_dir, "origin_by_group.tsv"), sep="\t")

    # --- function (eggNOG)
    e = pd.read_csv(os.path.join(a.outdir, "02_functional_annotation", "eggnog",
                                 "final_annotation_proteins_main.emapper.annotations"),
                    sep="\t", comment="#", header=None, dtype=str, low_memory=False)
    e = e.iloc[:, :13]
    e.columns = ["query", "seed", "evalue", "score", "OGs", "max_lvl", "COG", "Description",
                 "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway"]
    e["gene"] = e["query"].map(lambda q: CDS_SUFFIX.sub("", q))
    e = e.drop_duplicates("gene").set_index("gene")

    cog_rows, desc_rows, path_rows = [], [], []
    for g, d in groups.items():
        n = len(d)
        ed = e.reindex(d.index.intersection(e.index))
        cog = ed["COG"].dropna()
        cog = cog[~cog.isin(["-", ""])]
        cnt = pd.Series([c for s in cog for c in s]).value_counts()
        for c, k in cnt.items():
            cog_rows.append({"group": g, "COG": c, "category": COG.get(c, "?"), "n_genes": k,
                             "pct_of_group": 100 * k / n, "pct_of_annotated": 100 * k / len(cog)})
        desc = ed["Description"].dropna()
        desc = desc[~desc.isin(["-", ""])].value_counts().head(a.top)
        for x, k in desc.items():
            desc_rows.append({"group": g, "description": x, "n_genes": k, "pct_of_group": 100 * k / n})
        path = ed["KEGG_Pathway"].dropna()
        path = pd.Series([p for s in path[path != "-"] for p in s.split(",") if p.startswith("map")]).value_counts().head(a.top)
        for x, k in path.items():
            path_rows.append({"group": g, "kegg_pathway": x, "n_genes": k, "pct_of_group": 100 * k / n})
    cog_df = pd.DataFrame(cog_rows)
    cog_df.to_csv(os.path.join(a.out_dir, "function_cog.tsv"), sep="\t", index=False, float_format="%.2f")
    pd.DataFrame(desc_rows).to_csv(os.path.join(a.out_dir, "function_top_descriptions.tsv"), sep="\t",
                                   index=False, float_format="%.2f")
    pd.DataFrame(path_rows).to_csv(os.path.join(a.out_dir, "function_top_kegg.tsv"), sep="\t",
                                   index=False, float_format="%.2f")

    # --- markdown
    md = ["# Genes added by TITAN vs reference (Liftoff) genes\n",
          "Groups: " + ", ".join("%s = %d" % (g, len(d)) for g, d in groups.items()) + "\n",
          "## Support\n", support.round(1).set_index("group").T.to_string(),
          "\n## Expression (TPM)\n"]
    st = pd.DataFrame(stats)
    for m in ("max_tpm", "mean_tpm", "n_samples_expr"):
        sub = st[st.metric == m][["group", "n", "min", "q1", "median", "mean", "q3", "p90", "max"]].round(2)
        md += ["**%s**\n" % m, sub.to_string(index=False), ""]
    md += ["## Origin of the genes\n", by_origin.to_string(), "\n## eggNOG COG categories (% of group)\n"]
    if len(cog_df):
        pv = cog_df.pivot(index="category", columns="group", values="pct_of_group").fillna(0).round(1)
        md.append(pv.sort_values("new", ascending=False).to_string())
    open(os.path.join(a.out_dir, "summary.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
