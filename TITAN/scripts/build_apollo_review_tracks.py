#!/usr/bin/env python3
"""Write GFF3 tracks of questionable genes of a TITAN annotation, for manual review in Apollo / JBrowse.

NOT part of the pipeline: an analysis tool. Each flagged gene is copied with its mRNA / exon / CDS children,
and the gene line gets extra attributes (qc_flags, Note and one attribute per category) so the reason shows
in the Apollo feature panel. Categories:

  overlap      exons shared with another coding gene (classify_gene_overlaps.py --pairs-tsv); worst class
               per gene: cds_cds_same > cds_cds_opposite > cds_utr > utr_utr
  split        member of a candidate split-gene pair (estimate_split_genes.py, tiers high / likely)
  unsupported  evidence class other than expressed / conserved_or_functional (classify_unsupported_genes.py)
  monoexonic   single-exon gene with an unsupported_* tier (flag_low_confidence_monoexonic_genes.py)

Outputs in --out-dir: one <category>.gff3 per category, all_flagged.gff3 (union, sorted) and
flagged_genes.tsv (one row per flagged gene).
"""
import argparse
import csv
import os
from collections import defaultdict

OVERLAP_RANK = ["cds_cds_same", "cds_cds_opposite", "cds_utr", "utr_utr"]
UNSUPPORTED = {"no_signal", "liftoff_only", "junction_only", "te_protein", "te_overlap_unsupported"}
CATEGORIES = ["overlap", "split", "unsupported", "monoexonic"]


def attrs(raw):
    out = {}
    for item in raw.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def load_gff(path):
    """gene ID -> (gene line fields, [child lines]) and gene order."""
    genes, order, parent_gene = {}, [], {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9:
                continue
            a = attrs(f[8])
            if f[2] == "gene":
                genes[a["ID"]] = (f, [])
                order.append(a["ID"])
                continue
            parent = a.get("Parent", "").split(",")[0]
            g = parent if parent in genes else parent_gene.get(parent)
            if g is None:
                continue
            if "ID" in a:
                parent_gene[a["ID"]] = g
            genes[g][1].append(f)
    return genes, order


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff3", required=True)
    ap.add_argument("--overlap-pairs", required=True)
    ap.add_argument("--split-candidates", required=True)
    ap.add_argument("--gene-classes", required=True)
    ap.add_argument("--monoexonic", required=True)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    genes, order = load_gff(args.gff3)
    flags = defaultdict(dict)  # gene -> category -> attribute dict

    partners = defaultdict(lambda: defaultdict(list))
    for r in read_tsv(args.overlap_pairs):
        cls = r["class"] + "_" + r["strand"] if r["class"] == "cds_cds" else r["class"]
        partners[r["gene_a"]][cls].append(r["gene_b"])
        partners[r["gene_b"]][cls].append(r["gene_a"])
    for g, by_cls in partners.items():
        worst = min(by_cls, key=OVERLAP_RANK.index)
        flags[g]["overlap"] = {
            "overlap_class": worst,
            "overlap_partners": ",".join(sorted({p for ps in by_cls.values() for p in ps})),
        }

    for r in read_tsv(args.split_candidates):
        for g, other in ((r["gene_a"], r["gene_b"]), (r["gene_b"], r["gene_a"])):
            cur = flags[g].setdefault("split", {"split_tier": r["tier"], "split_partners": []})
            if r["tier"] == "high":
                cur["split_tier"] = "high"
            cur["split_partners"].append(other)
    for g in list(flags):
        if "split" in flags[g]:
            flags[g]["split"]["split_partners"] = ",".join(sorted(set(flags[g]["split"]["split_partners"])))

    for r in read_tsv(args.gene_classes):
        if r["evidence_class"] in UNSUPPORTED:
            flags[r["gene_id"]]["unsupported"] = {"evidence_class": r["evidence_class"]}

    for r in read_tsv(args.monoexonic):
        if r["tier"].startswith("unsupported"):
            flags[r["gene_id"]]["monoexonic"] = {"monoexonic_tier": r["tier"]}

    flags = {g: c for g, c in flags.items() if g in genes}
    os.makedirs(args.out_dir, exist_ok=True)
    handles = {c: open(os.path.join(args.out_dir, c + ".gff3"), "w") for c in CATEGORIES + ["all_flagged"]}
    for h in handles.values():
        h.write("##gff-version 3\n")
    counts = defaultdict(int)
    with open(os.path.join(args.out_dir, "flagged_genes.tsv"), "w") as tsv:
        tsv.write("gene_id\tchrom\tstart\tend\tstrand\tqc_flags\tdetails\n")
        for g in sorted(flags, key=lambda x: (genes[x][0][0], int(genes[x][0][3]))):
            gf, children = genes[g]
            extra = {}
            for c in CATEGORIES:
                extra.update(flags[g].get(c, {}))
            cats = [c for c in CATEGORIES if c in flags[g]]
            note = " | ".join(f"{k}: {v.replace(',', ' ')}" for k, v in extra.items())
            gline = gf[:8] + [gf[8].rstrip(";") + ";qc_flags=" + ",".join(cats) + ";Note=" + note + "".join(
                f";{k}={v}" for k, v in extra.items())]
            block = "\t".join(gline) + "\n" + "".join("\t".join(f) + "\n" for f in children) + "###\n"
            for c in cats + ["all_flagged"]:
                handles[c].write(block)
                counts[c] += 1
            tsv.write("\t".join([g, gf[0], gf[3], gf[4], gf[6], ",".join(cats), note]) + "\n")
    for h in handles.values():
        h.close()
    for c in CATEGORIES + ["all_flagged"]:
        print(f"{c}\t{counts[c]}")


if __name__ == "__main__":
    main()
