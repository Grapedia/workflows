#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
from tqdm import tqdm
from sklearn.mixture import GaussianMixture

# === ARGUMENT PARSING ===
# example : python3 scripts/filter_monoexons_from_GFF3.py --gff3 OUTDIR/augustus.hints.gff3 --output OUTDIR/augustus.filtered.gff3 --outdir OUTDIR/figures/
parser = argparse.ArgumentParser(description="Filter short mono-exonic transcripts from a GFF3 file.")
parser.add_argument("--gff3", required=True, help="Path to input GFF3 file")
parser.add_argument("--output", required=True, help="Path to output filtered GFF3 file")
parser.add_argument("--threshold", type=int, help="Minimum transcript length to retain (if not provided, automatically estimated)")
parser.add_argument("--outdir", type=str, help="Directory for output plots", default=".")
parser.add_argument("--quantile", type=float, default=10.0, help="Percentile to use for quantile threshold (default: 10)")
args = parser.parse_args()

# === PATH SETUP ===
gff3_path = Path(args.gff3)
output_filtered_gff = Path(args.output)
outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)  # Ensure output directory exists
output_plot = outdir / (output_filtered_gff.stem + ".png")  # Define path for plot image
print(f"Output directory: {outdir}")
print(f"Output plot: {output_plot}")

# === LOAD GFF3 ===
with open(gff3_path) as f:
    # Read all lines except comments and empty lines
    gff_lines = [line.strip() for line in f if not line.startswith("#") and line.strip()]

# === PARSE TO DATAFRAME ===
records = []
for line in gff_lines:
    fields = line.split('\t')
    if len(fields) != 9:
        continue  # Skip malformed lines
    chrom, source, feature, start, end, score, strand, phase, attributes = fields
    # Extract key=value pairs from the attributes column
    attr_dict = dict(re.findall(r'(\w+)=([^;]+)', attributes))
    records.append({
        "chrom": chrom, "source": source, "feature": feature,
        "start": int(start), "end": int(end),
        "score": score, "strand": strand, "phase": phase,
        "attributes": attributes,
        "ID": attr_dict.get("ID", None),
        "Parent": attr_dict.get("Parent", None),
        "gene_id": (
            attr_dict.get("Parent", None)  # for mRNA lines: Parent is the gene
            if feature == "mRNA"
            else attr_dict.get("Parent", None).split('.')[0]  # for exon/CDS: extract gene from transcript
            if attr_dict.get("Parent", None)
            else None
        )
    })

# Convert the list of dictionaries to a DataFrame
df = pd.DataFrame.from_records(records)
print("dataframe of the GFF3 file :")
print(df.head(10))

# === IDENTIFY mRNAs AND EXONS ===
mrnas = df[df["feature"] == "mRNA"].copy()
print("mRNAs dataframe created :")
print(mrnas.head(10))
exons = df[df["feature"] == "exon"].copy()
print("exons dataframe created :")
print(exons.head(10))

# Count how many exons each mRNA has
# exons["Parent"]: retrieves the Parent column of all exon → rows, i.e. the transcripts to which each exon belongs.
# .value_counts(): counts how many times each transcript appears = how many exons it has
# .to_dict(): converts the result into a dictionary
# Example:
# {
#  "g1.t1": 3,
#  "g2.t1": 1,
#  "g3.t1": 2
# }
exon_counts = exons["Parent"].value_counts().to_dict()
# mrnas["ID"]: these are the identifiers of the mRNAs (e.g. "g1.t1")
# .map(exon_counts): we associate the number of exons with each mRNA, using the dictionary
# .fillna(0): if an mRNA has no associated exon (rare but possible), we set this to 0
# .astype(int): to ensure that the values are integers
mrnas["exon_count"] = mrnas["ID"].map(exon_counts).fillna(0).astype(int)

# Compute transcript length
mrnas["length"] = mrnas["end"] - mrnas["start"] + 1

# Identify genes that have only one transcript
gene_transcript_counts = mrnas["gene_id"].value_counts()
unique_transcript_genes = gene_transcript_counts[gene_transcript_counts == 1].index

print("mRNAs dataframe with exon_count and length added :")
print(mrnas.head(10))

# Keep only mRNAs that are mono-exonic and from single-transcript genes
monoexonic_mrnas = mrnas[
    (mrnas["exon_count"] == 1) &
    (mrnas["gene_id"].isin(unique_transcript_genes))
].copy()

# Exit early if no mono-exonic transcripts are found
if monoexonic_mrnas.empty:
    print("No mono-exonic transcripts found in the input file. No filtering applied.")
    Path(output_filtered_gff).write_text("".join(Path(gff3_path).read_text().splitlines()) + "")
    exit(0)

# === AUTOMATIC THRESHOLD ESTIMATION IF NOT PROVIDED ===
if args.threshold is not None:
    final_threshold = args.threshold  # Use user-defined threshold
else:
    lengths = monoexonic_mrnas["length"].values
    lengths_sorted = np.sort(lengths)
    x = np.arange(len(lengths_sorted))

    # Method 1: Quantile (e.g. 10th percentile)
    quantile_threshold = np.percentile(lengths, args.quantile)

    # Method 2: GMM
    try:
        gmm = GaussianMixture(n_components=2, random_state=0).fit(lengths.reshape(-1, 1))
        means = gmm.means_.flatten()
        gmm_threshold = min(means)
    except:
        gmm_threshold = np.inf

    # Final threshold = minimum of all methods (most conservative)
    final_threshold = int(min(quantile_threshold, gmm_threshold))


print(f"Method 1 : Quantile (10th percentile), threshold found : {quantile_threshold} bp")
print(f"Method 2 : GMM method, threshold found : {gmm_threshold} bp")

print(f"Chosen threshold (minimum of all methods, most conservative): {final_threshold} bp")

# === PLOT LENGTH DISTRIBUTION ===
plt.figure(figsize=(10, 6))
plt.hist(monoexonic_mrnas["length"], bins=50, color='skyblue', edgecolor='black')
plt.axvline(final_threshold, color='red', linestyle='dashed', label=f"Final threshold = {final_threshold} bp")
plt.axvline(quantile_threshold, color='green', linestyle='--', label=f"{args.quantile}th percentile")
plt.axvline(gmm_threshold, color='brown', linestyle='--', label='GMM')
plt.title("Length distribution of unique mono-exonic transcripts")
plt.xlabel("Transcript length (bp)")
plt.ylabel("Number of transcripts")
plt.legend()
plt.tight_layout()
plt.savefig(output_plot)

# === FILTER GFF3 ===
# Collect IDs of transcripts and genes to be removed based on threshold
to_remove_ids = monoexonic_mrnas[monoexonic_mrnas["length"] < final_threshold]["ID"].tolist()
to_remove_gene_ids = monoexonic_mrnas[monoexonic_mrnas["length"] < final_threshold]["Parent"].tolist()

def keep_line(line):
    if line.startswith("#"):
        return True
    fields = line.strip().split('\t')
    if len(fields) != 9:
        return True
    attr_field = fields[8]
    attr_dict = dict(re.findall(r'(\w+)=([^;]+)', attr_field))
    feature_id = attr_dict.get("ID", "")
    parent_id = attr_dict.get("Parent", "")
    
    # Remove lines if the feature is a transcript or exon from a removed transcript
    if feature_id in to_remove_ids or parent_id in to_remove_ids:
        return False
    
    # Remove gene feature if its ID is in the list of removed gene IDs
    if fields[2] == "gene" and feature_id in to_remove_gene_ids:
        return False

    return True

# Apply filtering to original GFF3 content
with open(gff3_path) as fin, open(output_filtered_gff, "w") as fout:
    for line in tqdm(fin, desc="Filtering GFF3"):
        if keep_line(line):
            fout.write(line)

# Identify elements removed
removed_transcripts = monoexonic_mrnas[monoexonic_mrnas["length"] < final_threshold]
removed_transcript_ids = set(removed_transcripts["ID"])
removed_gene_ids = set(removed_transcripts["gene_id"])

# From full GFF dataframe
removed_exons = df[
    (df["feature"] == "exon") &
    (df["Parent"].isin(removed_transcript_ids))
]

# Summary printout
print("\n=== Overview ===")
print(f"Total genes in GFF3                        : {df[df['feature'] == 'gene']['ID'].nunique()}")
print(f"Total transcripts (mRNAs)                 : {len(mrnas)}")
print(f"Total of transcripts that are mono-exonic and from single-transcript genes : {len(monoexonic_mrnas)}")

print(f"Summary of elements removed below {final_threshold} bp:")
print(f"- Genes removed     : {len(removed_gene_ids)}")
print(f"- Transcripts removed: {len(removed_transcript_ids)}")
print(f"- Exons removed     : {len(removed_exons)}")

print(f"Conclusion : of the {len(monoexonic_mrnas)} monoexonic transcripts, {len(removed_transcript_ids)} are deleted from gff3 because they are too short")

print(f"Filtered GFF3 saved to: {output_filtered_gff}")
print(f"Plot saved to: {output_plot}")
