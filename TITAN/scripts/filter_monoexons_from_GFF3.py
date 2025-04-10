#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np
from sklearn.cluster import KMeans
from tqdm import tqdm
from sklearn.mixture import GaussianMixture
from skimage.filters import threshold_otsu

# === ARGUMENT PARSING ===
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
print("exons dataframe created :")
exons = df[df["feature"] == "exon"].copy()
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

    # Method 21: KMeans (bimodal clustering)
    kmeans = KMeans(n_clusters=2, random_state=0).fit(lengths.reshape(-1, 1))
    cluster_labels = kmeans.labels_
    cluster_0_max = lengths[cluster_labels == 0].max()
    cluster_1_max = lengths[cluster_labels == 1].max()
    kmeans_threshold = min(cluster_0_max, cluster_1_max)

    # Method 2: Quantile (e.g. 10th percentile)
    quantile_threshold = np.percentile(lengths, args.quantile)

    # Method 3: Otsu's method
    try:
        otsu_threshold = threshold_otsu(lengths)
    except:
        otsu_threshold = np.inf

    # Method 4: GMM
    try:
        gmm = GaussianMixture(n_components=2, random_state=0).fit(lengths.reshape(-1, 1))
        means = gmm.means_.flatten()
        gmm_threshold = min(means)
    except:
        gmm_threshold = np.inf

    # Final threshold = minimum of all methods (most conservative)
    final_threshold = int(min(kmeans_threshold, quantile_threshold, otsu_threshold, gmm_threshold))


print(f"Method 1 : KMeans (bimodal clustering), threshold found : {kmeans_threshold} bp")
print(f"Method 2 : Quantile (10th percentile), threshold found : {quantile_threshold} bp")
print(f"Method 3 : Otsu's method, threshold found : {otsu_threshold} bp")
print(f"Method 4 : GMM method, threshold found : {gmm_threshold} bp")

print(f"Chosen threshold (minimum of all methods, most conservative): {final_threshold} bp")

# === PLOT LENGTH DISTRIBUTION ===
plt.figure(figsize=(10, 6))
plt.hist(monoexonic_mrnas["length"], bins=50, color='skyblue', edgecolor='black')
plt.axvline(final_threshold, color='red', linestyle='dashed', label=f"Final threshold = {final_threshold} bp")
plt.axvline(kmeans_threshold, color='purple', linestyle='--', label='KMeans')
plt.axvline(quantile_threshold, color='green', linestyle='--', label=f"{args.quantile}th percentile")
plt.axvline(otsu_threshold, color='blue', linestyle='--', label='Otsu')
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
        return True  # Keep comment lines
    fields = line.strip().split('\t')
    if len(fields) != 9:
        return True  # Keep malformed lines just in case
    attr_field = fields[8]
    attr_dict = dict(re.findall(r'(\w+)=([^;]+)', attr_field))
    feature_id = attr_dict.get("ID", "")
    parent_id = attr_dict.get("Parent", "")
    # Discard lines with IDs matching those to remove
    if feature_id in to_remove_ids or parent_id in to_remove_ids:
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
print(f"Summary of elements removed below {final_threshold} bp:")
print(f"- Genes removed     : {len(removed_gene_ids)}")
print(f"- Transcripts removed: {len(removed_transcript_ids)}")
print(f"- Exons removed     : {len(removed_exons)}")

print(f"Filtered GFF3 saved to: {output_filtered_gff}")
print(f"Plot saved to: {output_plot}")
