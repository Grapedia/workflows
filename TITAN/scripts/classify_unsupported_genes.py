#!/usr/bin/env python3
"""Classify genes with no RNA-seq expression by what *else* supports them.

Reads gene_evidence_table.tsv (build_gene_evidence_table.py + add_intron_support.py +
add_homology_features.py) and assigns each protein-coding gene one `evidence_class`:

  expressed                TPM >= 0.5 in >= 1 of the RNA-seq libraries (reference set)
  te_protein               not expressed AND >= 50 % of the CDS lies in an EDTA TE annotation AND a
                           functional hit describes a transposable-element protein (RT, gag-pol,
                           integrase, transposase ...): a real TE-encoded protein, not a host gene
  te_overlap_unsupported   not expressed, >= 50 % of the CDS in an EDTA TE annotation, no TE-protein
                           hit, and no function / OMAMER / non-Vitis homolog: a gene call on TE sequence
  conserved_or_functional  not expressed, not one of the two above, but has a functional hit
                           (Diamond2GO / eggNOG / InterProScan / MapMan), an OMAMER hog placement,
                           or a >=35 % identity / >=50 % coverage homolog in a non-Vitis eudicot
  liftoff_only             none of the above, but carried over from the v4.3 annotation
  junction_only            none of the above, but >=1 intron is reproduced by an RNA-seq assembly
  no_signal                nothing at all: no expression, no junction, no function, no orthology,
                           no non-Vitis homolog, not a v4.3 gene, not a TE

Independent checks (never used to define the classes) are added as columns:
  in_v51 (>=50 % of the CDS overlaps a PN40024 5.1 CDS, same strand), busco (BUSCO ortholog).

For `no_signal`/`te_overlap_unsupported`/`liftoff_only`, artefact signatures explain *why* the model is
suspect (see docs/user/unsupported_genes.md).  Outputs (in --out-dir):
  gene_classes.tsv            one row per gene: class, signatures, key evidence
  class_summary.tsv           counts per class and per origin
  signature_summary.tsv       signature prevalence per class
  set_aside_candidates.txt    gene IDs proposed to be set aside (see --policy)
  summary.json
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

MIN_HOMOL_PIDENT = 35.0
MIN_HOMOL_QCOV = 0.5
TE_CDS_FRAC = 0.5
SHORT_PROT = 100


def load(table, v51_gff_overlap, busco_table):
    d = pd.read_csv(table, sep="\t", index_col=0)
    if v51_gff_overlap and os.path.exists(v51_gff_overlap):
        ov = pd.read_csv(v51_gff_overlap, sep="\t", names=["g", "g51", "bp"]).groupby("g").bp.max()
        d["v51_cds_frac"] = (ov.reindex(d.index).fillna(0) / d["cds_len"].clip(lower=1)).clip(upper=1)
        d["in_v51"] = d["v51_cds_frac"] >= 0.5
    if busco_table and os.path.exists(busco_table):
        b = pd.read_csv(busco_table, sep="\t", comment="#", header=None)
        b = b[b[1].isin(["Complete", "Duplicated", "Fragmented"])]
        ids = set(b[2].str.replace(r"_t\d+_CDS\d+\.prot$", "", regex=True))
        d["busco"] = d.index.isin(ids)
    return d


def classify(d):
    d = d[d["biotype"] == "mRNA"].copy()
    d["has_function"] = d["functional_any"] | d["mapman_any"]
    d["has_junction"] = (d["n_introns"] > 0) & (d["frac_introns_any"] > 0)
    d["has_nonvitis_homolog"] = ((d["ref_nonvitis_qcov"] >= MIN_HOMOL_QCOV) &
                                 (d["ref_nonvitis_pident"] >= MIN_HOMOL_PIDENT))
    d["te_overlap"] = d["te_frac_cds"] >= TE_CDS_FRAC
    # te_kind is defined for every gene (expressed or not): the TE flag does not depend on expression
    d["te_kind"] = np.select(
        [d["te_overlap"] & d["te_domain_hit"], d["te_overlap"], d["te_domain_hit"]],
        ["te_protein", "te_overlap_only", "te_domain_only"], "none")
    supported = d["has_function"] | d["omamer_hit"] | d["has_nonvitis_homolog"]
    d["evidence_class"] = np.select(
        [d["expressed"],
         d["te_kind"] == "te_protein",
         (d["te_kind"] == "te_overlap_only") & ~supported,
         supported,
         d["liftoff_carried"],
         d["has_junction"]],
        ["expressed", "te_protein", "te_overlap_unsupported", "conserved_or_functional",
         "liftoff_only", "junction_only"],
        "no_signal")
    return d


def signatures(d):
    """Artefact signatures (boolean columns).  Descriptive only: not used to build classes."""
    s = pd.DataFrame(index=d.index)
    s["sig_short_orf"] = d["prot_len"] < SHORT_PROT
    s["sig_single_exon"] = d["single_exon"]
    s["sig_raw_abinitio_not_tsebra"] = d["origin"].isin(["ab_initio_braker"]) & ~d["in_tsebra"]
    s["sig_helixer_only"] = d["origin"] == "ab_initio_helixer"
    s["sig_te_overlap"] = d["te_frac_cds"] >= TE_CDS_FRAC
    s["sig_te_domain"] = d["te_domain_hit"]
    s["sig_fragment_of_paralog"] = ((d["self_qcov"] < 0.9) & (d["self_scov"] >= 0.5) &
                                    (d["self_pident"] >= 50))
    s["sig_near_supported_same_strand"] = d["near_supported_same_strand"].fillna(False).astype(bool)
    s["sig_antisense_or_nested"] = (d["ovl_antisense"] > 0) | (d["nested_in_other"] > 0)
    s["sig_cds_overlap_same_strand"] = d["ovl_cds_same_strand"] > 0
    s["sig_low_complexity"] = (d["prot_entropy"] < 3.5) | (d["prot_top_aa_frac"] > 0.2)
    s["sig_no_orthology"] = ~d["omamer_hit"]
    s["sig_chr00"] = d["chrom"] == "chr00"
    s["sig_identical_copy_elsewhere"] = (d["self_pident"] >= 95) & (d["self_qcov"] >= 0.8)
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--v51-overlap", help="TSV gene_final<TAB>gene_v51<TAB>overlapping_cds_bp")
    ap.add_argument("--busco-table")
    ap.add_argument("--salmon-numreads", help="TSV from salmon_gene_numreads.py")
    ap.add_argument("--policy", default="strict", choices=["strict", "broad"],
                    help="strict: no_signal + te_overlap_unsupported; broad: + liftoff_only + junction_only")
    a = ap.parse_args()

    d = load(a.table, a.v51_overlap, a.busco_table)
    if a.salmon_numreads:
        d = d.join(pd.read_csv(a.salmon_numreads, sep="\t", index_col=0))
    d = classify(d)
    sig = signatures(d)
    d = d.join(sig)
    os.makedirs(a.out_dir, exist_ok=True)

    keep = ["chrom", "start", "end", "strand", "origin", "mikado_alias", "n_exons", "cds_len", "prot_len",
            "max_tpm", "n_samples_expr", "n_samples_any", "evidence_class", "has_function", "omamer_hit",
            "has_nonvitis_homolog", "has_junction", "liftoff_carried", "in_tsebra", "te_frac_cds",
            "te_domain_hit", "te_kind", "self_pident", "self_qcov"] + list(sig.columns)
    keep += [c for c in ("in_v51", "busco", "salmon_reads_total", "salmon_n_samples_ge10reads") if c in d.columns]
    d[[c for c in keep if c in d.columns]].to_csv(os.path.join(a.out_dir, "gene_classes.tsv"), sep="\t")

    order = ["expressed", "conserved_or_functional", "te_protein", "te_overlap_unsupported",
             "liftoff_only", "junction_only", "no_signal"]
    cs = d["evidence_class"].value_counts().reindex(order).to_frame("n_genes")
    cs["pct_of_protein_coding"] = (100 * cs["n_genes"] / len(d)).round(2)
    for col in ("in_v51", "busco"):
        if col in d.columns:
            cs[col] = d.groupby("evidence_class")[col].sum().reindex(order)
    if "in_v51" in d.columns:
        cs["pct_in_v51"] = (100 * d.groupby("evidence_class")["in_v51"].mean().reindex(order)).round(1)
    orig = pd.crosstab(d["evidence_class"], d["origin"]).reindex(order)
    pd.concat([cs, orig], axis=1).to_csv(os.path.join(a.out_dir, "class_summary.tsv"), sep="\t")

    if "salmon_reads_total" in d.columns:
        rd = d.groupby("evidence_class")["salmon_reads_total"].agg(
            median_reads="median", pct_zero_reads=lambda x: 100 * (x == 0).mean(),
            pct_lt10_reads=lambda x: 100 * (x < 10).mean(), pct_ge100_reads=lambda x: 100 * (x >= 100).mean()
        ).reindex(order).round(1)
        rd.to_csv(os.path.join(a.out_dir, "salmon_reads_by_class.tsv"), sep="\t")
    sp = d.groupby("evidence_class")[list(sig.columns)].mean().reindex(order).T.round(3)
    sp.to_csv(os.path.join(a.out_dir, "signature_summary.tsv"), sep="\t")

    classes = (["no_signal", "te_overlap_unsupported"] if a.policy == "strict"
               else ["no_signal", "te_overlap_unsupported", "liftoff_only", "junction_only"])
    aside = d[d["evidence_class"].isin(classes)]
    with open(os.path.join(a.out_dir, "set_aside_candidates.txt"), "w") as fh:
        fh.write("\n".join(aside.index) + "\n")
    summ = {"protein_coding_genes": int(len(d)),
            "not_expressed": int((~d["expressed"]).sum()),
            "class_counts": {k: int(v) for k, v in cs["n_genes"].items()},
            "policy": a.policy, "set_aside_candidates": int(len(aside)),
            "thresholds": {"min_tpm": 0.5, "homolog_pident": MIN_HOMOL_PIDENT, "homolog_qcov": MIN_HOMOL_QCOV,
                           "te_cds_frac": TE_CDS_FRAC, "short_orf_aa": SHORT_PROT}}
    with open(os.path.join(a.out_dir, "summary.json"), "w") as fh:
        json.dump(summ, fh, indent=2)
    print(cs.to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
