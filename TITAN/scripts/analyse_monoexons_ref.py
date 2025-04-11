#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt

# example : python3 scripts/analyse_monoexons_ref.py --gff3 data/annotations/v4_3_just_ref.gff3 --outplot v4_3_just_ref_monoexonic_lengths.png
parser = argparse.ArgumentParser(description="Analyze transcript lengths for mono-exonic transcripts (reference GFF3)")
parser.add_argument("--gff3", required=True, help="Reference GFF3 file")
parser.add_argument("--outplot", default="monoexonic_lengths_ref.png", help="Output filename for the histogram")
args = parser.parse_args()

# === Read GFF3 ===
with open(args.gff3) as f:
    gff_lines = [line.strip() for line in f if not line.startswith("#") and line.strip()]

records = []
for line in gff_lines:
    fields = line.split('\t')
    if len(fields) != 9:
        continue
    chrom, source, feature, start, end, score, strand, phase, attributes = fields
    attr_dict = dict(re.findall(r'(\w+)=([^;]+)', attributes))
    if feature == "gene":
        gene_id = attr_dict.get("ID", None)
    elif feature == "mRNA":
        gene_id = attr_dict.get("Parent", None)
    elif feature == "exon":
        gene_id = None  # optional: you could map this later if needed
    else:
        gene_id = None
    records.append({
        "chrom": chrom, "source": source, "feature": feature,
        "start": int(start), "end": int(end),
        "strand": strand, "attributes": attributes,
        "ID": attr_dict.get("ID", None),
        "Parent": attr_dict.get("Parent", None),
        "gene_id": gene_id
    })

df = pd.DataFrame(records)

# Count total number of genes (before filtering)
total_genes = df[df["feature"] == "gene"]["ID"].nunique()

df = df[df["feature"].isin(["mRNA", "exon"])]
print("dataframe of the GFF3 file :")
print(df.head(10))

# Extract mRNAs and exons
mrnas = df[df["feature"] == "mRNA"].copy()
print("mRNAs dataframe created :")
print(mrnas.head(10))
exons = df[df["feature"] == "exon"].copy()
print("exons dataframe created :")
print(exons.head(10))

# Count number of exons per transcript
exon_counts = exons["Parent"].value_counts().to_dict()
mrnas["exon_count"] = mrnas["ID"].map(exon_counts).fillna(0).astype(int)

# Compute transcript length
mrnas["length"] = mrnas["end"] - mrnas["start"] + 1

# Keep only mono-exonic transcripts from single-transcript genes
gene_transcript_counts = mrnas["gene_id"].value_counts()
unique_transcript_genes = gene_transcript_counts[gene_transcript_counts == 1].index
monoexonic_mrnas = mrnas[
    (mrnas["exon_count"] == 1) &
    (mrnas["gene_id"].isin(unique_transcript_genes))
].copy()

# Analyze length distribution
lengths = monoexonic_mrnas["length"]
print("==== Statistics for mono-exonic transcripts ====")
print(f"Total genes                   : {total_genes}")
print(f"Total mono-exonic transcripts : {len(lengths)}")
print(f"Minimum length                : {lengths.min()} bp")
print(f"Maximum length                : {lengths.max()} bp")
print(f"Mean length                   : {lengths.mean():.1f} bp")
print(f"Median length                 : {lengths.median()} bp")
print(f"10th percentile               : {np.percentile(lengths, 10):.1f} bp")
print(f"5th percentile                : {np.percentile(lengths, 5):.1f} bp")
print(f"1st percentile                : {np.percentile(lengths, 1):.1f} bp")

# Plot length distribution
plt.figure(figsize=(10, 5))
plt.hist(lengths, bins=50, color="skyblue", edgecolor="black")
plt.axvline(np.percentile(lengths, 10), color="green", linestyle="--", label="10th percentile")
plt.axvline(lengths.min(), color="red", linestyle="--", label="Minimum observed")
plt.xlabel("Transcript length (bp)")
plt.ylabel("Number of transcripts")
plt.title("Length distribution of mono-exonic transcripts (reference)")
plt.legend()
plt.tight_layout()
plt.savefig(args.outplot)
print(f"Histogram saved to: {args.outplot}")
