#!/bin/bash
# Reproduce the "unsupported genes" analysis of a finished TITAN run (see docs/user/unsupported_genes.md).
#
#   scripts/run_unsupported_genes_analysis.sh data/titan_prod_out data/unsupported_genes_analysis [v5.1.gff3]
#
# Steps 1-3 need python3 + pandas + bedtools.  Step 2 (DIAMOND) needs the diamond2go image and the
# protein_data/ reference FASTAs; run it on a compute node (see diamond/run_diamond.sh).
set -euo pipefail
OUT=${1:?TITAN output dir}; WORK=${2:?analysis dir}; V51=${3:-}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$WORK"
TABLE=$WORK/gene_evidence_table.tsv

echo "[1/4] gene evidence table";           python3 "$HERE/build_gene_evidence_table.py" --outdir "$OUT" --out "$TABLE"
python3 "$HERE/add_intron_support.py" --outdir "$OUT" --table "$TABLE"

echo "[2/4] homology (DIAMOND) -- expects $WORK/diamond/{self,ref_vs,ref_eu}.tsv"
for f in self ref_vs ref_eu; do [[ -s $WORK/diamond/$f.tsv ]] || { echo "missing $WORK/diamond/$f.tsv: run diamond/run_diamond.sh first" >&2; exit 1; }; done
python3 "$HERE/add_homology_features.py" --outdir "$OUT" --table "$TABLE" --diamond-dir "$WORK/diamond"

V51ARG=""
if [[ -n $V51 ]]; then
  echo "[3/4] overlap with an independent reference annotation ($V51)"
  bed() { awk -F'\t' '$3=="CDS"{split($9,a,"Parent=");sub(/;.*/,"",a[2]);sub(/_t[0-9]+$/,"",a[2]);print $1"\t"$4-1"\t"$5"\t"a[2]"\t.\t"$7}' "$1" | sort -k1,1 -k2,2n; }
  bed "$OUT/01_final_annotation/primary/final_annotation.gff3" > "$WORK/fin_cds.bed"
  bed "$V51" > "$WORK/v51_cds.bed"
  bedtools intersect -s -wo -a "$WORK/fin_cds.bed" -b "$WORK/v51_cds.bed" 2>/dev/null \
    | awk -F'\t' '{k=$4"\t"$10; d[k]+=$13} END{for(k in d) print k"\t"d[k]}' > "$WORK/fin_vs_v51.tsv"
  V51ARG="--v51-overlap $WORK/fin_vs_v51.tsv"
fi

echo "[4/4] Salmon read totals + classification"
python3 "$HERE/salmon_gene_numreads.py" --quants-dir "$OUT/01_final_annotation/quality_report/expression_validation/quants" --out "$WORK/salmon_numreads.tsv"
python3 "$HERE/classify_unsupported_genes.py" --table "$TABLE" --out-dir "$WORK/classes" $V51ARG \
  --salmon-numreads "$WORK/salmon_numreads.tsv" --busco-table "$OUT/01_final_annotation/quality_report/busco/busco_full_table.tsv"
