#!/usr/bin/env python3
"""List and remove the TE-derived genes of a TITAN GFF3 that nothing else supports.

Post-hoc analysis tool (same spirit as flag_unsupported_genes.py): it works from the evidence table of
scripts/run_unsupported_genes_analysis.sh and never modifies its inputs or renames genes.  The pipeline itself
runs scripts/filter_te_ncrna_final_annotation.py (modules/primary_noTE.nf), its in-pipeline version.

A gene is removed when ALL of these hold
  1. >= --te-frac (default 0.80) of its CDS lies in an EDTA TE            (te_frac_cds)
  2. it was not carried over from the Liftoff reference annotation        (liftoff_carried)
  3. its CDS does not overlap the PN40024 5.1 reference annotation        (in_v51, gene_classes.tsv)
  4. no non-Vitis eudicot homolog (>= 35 % identity, >= 50 % query cover)  (ref_nonvitis_*)
  5. no functional hit other than a TE protein (Diamond2GO GO, eggNOG description/GO/KEGG,
     InterProScan informative analysis)                                    (functional_strict, te_domain_hit)
Conditions 2-5 protect real host genes that merely contain a TE insertion.  Run it after
scripts/run_unsupported_genes_analysis.sh, with InterProScan and MapMan results available: without them
condition 5 is too permissive.

Inputs: gene_evidence_table.tsv and classes/gene_classes.tsv of scripts/run_unsupported_genes_analysis.sh.

Usage
  filter_te_genes.py --gff3 final_annotation.gff3 \\
      --table data/unsupported_genes_analysis/gene_evidence_table.tsv \\
      --classes data/unsupported_genes_analysis/classes/gene_classes.tsv \\
      --out final_annotation.noTE.gff3 --removed-list te_genes_removed.tsv \\
      [--proteins final_annotation_proteins_main.fasta --proteins-out proteins_main.noTE.fasta] [--summary-only]

Python >= 3.6, pandas required.
"""
import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flag_unsupported_genes import dropped_feature_ids, filter_fasta, gff_records, parse_attrs  # noqa: E402

MIN_HOMOL_PIDENT = 35.0   # same thresholds as classify_unsupported_genes.py
MIN_HOMOL_QCOV = 0.5


def as_bool(s, index):
    return s.reindex(index).astype(str).eq("True")


def select(table, classes, te_frac):
    t = pd.read_csv(table, sep="\t", index_col="gene_id")
    c = pd.read_csv(classes, sep="\t", index_col="gene_id")
    te = t["te_frac_cds"].fillna(0) >= te_frac
    liftoff = as_bool(t["liftoff_carried"], t.index)
    in_v51 = as_bool(c["in_v51"], t.index)
    nonvitis = (t["ref_nonvitis_qcov"].fillna(0) >= MIN_HOMOL_QCOV) & (t["ref_nonvitis_pident"].fillna(0) >= MIN_HOMOL_PIDENT)
    host_function = as_bool(t["functional_strict"], t.index) & ~as_bool(t["te_domain_hit"], t.index)
    t["remove"] = te & ~liftoff & ~in_v51 & ~nonvitis & ~host_function
    return t


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff3", required=True, help="final_annotation.gff3")
    ap.add_argument("--table", required=True, help="gene_evidence_table.tsv")
    ap.add_argument("--classes", required=True, help="classes/gene_classes.tsv")
    ap.add_argument("--te-frac", type=float, default=0.80, help="minimum fraction of CDS in a TE (default 0.80)")
    ap.add_argument("--out", help="filtered GFF3")
    ap.add_argument("--removed-list", help="TSV of the removed genes with their evidence")
    ap.add_argument("--proteins", help="protein FASTA with CDS-record headers to filter too")
    ap.add_argument("--proteins-out")
    ap.add_argument("--summary-only", action="store_true")
    a = ap.parse_args()

    t = select(a.table, a.classes, a.te_frac)
    drop = set(t.index[t["remove"]])
    r = t[t["remove"]]
    print("genes in table: %d; TE-derived genes to remove (CDS in TE >= %.0f %%): %d (%.1f %%)"
          % (len(t), 100 * a.te_frac, len(drop), 100 * len(drop) / len(t)))
    print("  expressed: %.1f %%; median max TPM %.2f; >= 3 samples: %d; max TPM >= 10: %d"
          % (100 * r["expressed"].mean(), r["max_tpm"].median(), (r["n_samples_expr"] >= 3).sum(), (r["max_tpm"] >= 10).sum()))
    print("  by origin: " + ", ".join("%s %d" % kv for kv in r["origin"].value_counts().items()))
    if a.removed_list:
        cols = ["origin", "chrom", "start", "end", "strand", "te_frac_cds", "max_tpm", "n_samples_expr",
                "functional_strict", "te_domain_hit", "ref_nonvitis_pident", "ref_nonvitis_qcov"]
        r[cols].to_csv(a.removed_list, sep="\t")
    if a.summary_only or not a.out:
        return

    dropped = dropped_feature_ids(a.gff3, drop)
    kept = removed_lines = 0
    with open(a.out, "w") as out:
        for line, f in gff_records(a.gff3):
            if f is None:
                out.write(line)
                continue
            attrs = parse_attrs(f[8])
            parents = attrs.get("Parent", "").split(",") if "Parent" in attrs else []
            if attrs.get("ID") in dropped or (parents and all(p in dropped for p in parents)):
                removed_lines += 1
                continue
            kept += f[2] == "gene"
            out.write(line)
    print("GFF3: %d genes kept, %d lines removed -> %s" % (kept, removed_lines, a.out))
    if a.proteins:
        filter_fasta(a.proteins, a.proteins_out, drop)
        print("proteins -> %s" % a.proteins_out)


if __name__ == "__main__":
    main()
