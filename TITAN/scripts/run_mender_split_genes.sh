#!/usr/bin/env bash
# Detect and merge split gene models with Mender (https://github.com/simrbase/mender).
#
# Mender aligns the annotation proteins against a reference proteome (DIAMOND blastp) and calls
# neighbouring genes whose alignments tile end to end on the same reference protein; Iso-Seq reads
# spanning the fragments are extra evidence. Merges are applied to the GFF3, then each merged protein
# is translated and its junction scored against reference proteins (PASS / REVIEW / FAIL); merges
# that do not pass go back to their source genes. This wrapper:
#   1. renames the TITAN protein headers (<transcript>_CDSn.prot) to transcript IDs for Mender, and
#      gives Mender only the genes with an mRNA (it rejects ncRNA transcripts), with exons shared by
#      several isoforms (Parent=t1,t2) written once per transcript (Mender rejects multi-parent features),
#   2. builds the reference from the protein evidence FASTAs, minus records of --exclude-os
#      (OS=<genus> in UniProt headers): a reference built from the same genus would repeat its errors,
#   3. turns long-read BAMs into one mRNA feature per primary alignment (MAPQ >= 20),
#   4. writes the Mender config and runs steps 1-8 and 10 (AGAT step 9 is skipped: TITAN runs AGAT),
#   5. writes the merged annotation in TITAN format: final_annotation.gff3 = the input GFF3 minus the
#      source genes of the validated merges plus the merged genes, in genome order (every other gene is
#      kept exactly as it was); proteins_all / proteins_main likewise (merged transcripts translated
#      from the genome),
#      plus mender_merges.tsv (one row per merged gene) and mender_summary.json.
#
# Usage: run_mender_split_genes.sh --gff3 F --proteins-all F --proteins-main F --genome F
#          --reference "F1 F2 ..." [--long-read-bams "B1 B2 ..."] [--threads N] [--evalue E]
#          [--exclude-os GENUS] [--flags L] [--skip-flags L] [--include-review yes|no]
#          [--gene-template T] [--trans-template T]
set -euo pipefail

gff3=""; proteins_all=""; proteins_main=""; genome=""; reference=""; bams=""
threads=4; evalue="1e-20"; exclude_os="Vitis"; flags="STRONG,CLEAN"
skip_flags="SKIPPED_GENE,LOW_COV,SINGLE_HIT,TRANSITIVE_JOIN,LARGE_SPAN_EXTREME"
include_review="no"; gene_template="VitviMERGEg[GCOUNT:5]"; trans_template="VitviMERGEg[GCOUNT:5]_t[TCOUNT:3]"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --gff3) gff3="$2" ;;
    --proteins-all) proteins_all="$2" ;;
    --proteins-main) proteins_main="$2" ;;
    --genome) genome="$2" ;;
    --reference) reference="$2" ;;
    --long-read-bams) bams="$2" ;;
    --threads) threads="$2" ;;
    --evalue) evalue="$2" ;;
    --exclude-os) exclude_os="$2" ;;
    --flags) flags="$2" ;;
    --skip-flags) skip_flags="$2" ;;
    --include-review) include_review="$2" ;;
    --gene-template) gene_template="$2" ;;
    --trans-template) trans_template="$2" ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
  shift 2
done
for v in gff3 proteins_all proteins_main genome reference; do
  [[ -n "${!v}" ]] || { echo "Missing --${v//_/-}" >&2; exit 2; }
done
mender_dir="${MENDER_DIR:-/opt/mender}"

# 1. Mender input: coding genes only (it rejects ncRNA transcripts), one line per parent for features
#    shared by several isoforms (it rejects Parent=t1,t2), and query proteins with transcript-ID headers
python3 - "$gff3" <<'PY'
import re, sys
blocks, cur = [], None
for line in open(sys.argv[1]):
    f = line.rstrip("\n").split("\t")
    if len(f) != 9:
        continue
    if f[2] == "gene":
        cur = [[line], False]
        blocks.append(cur)
    elif cur is not None:
        parent = re.search(r"(?:^|;)Parent=([^;]+)", f[8])
        parents = parent.group(1).split(",") if parent else []
        if len(parents) > 1:
            for p in parents:
                attrs = f[8].replace(parent.group(0), parent.group(0).replace(parent.group(1), p), 1)
                attrs = re.sub(r"(^|;)ID=([^;]+)", lambda m: f"{m.group(1)}ID={m.group(2)}-{p}", attrs, count=1)
                cur[0].append("\t".join(f[:8] + [attrs]) + "\n")
        else:
            cur[0].append(line)
        if f[2] == "mRNA":
            cur[1] = True
with open("coding_genes.gff3", "w") as c:
    c.write("##gff-version 3\n")
    for lines, coding in blocks:
        if coding:
            c.writelines(lines + ["###\n"])
PY
sed -E '/^>/ s/^>([^ ]+)_CDS[0-9]+\.prot.*/>\1/' "$proteins_main" > query_proteins.fasta
awk -F'\t' '$3=="mRNA"{match($9,/ID=[^;]+/); print substr($9,RSTART+3,RLENGTH-3)}' coding_genes.gff3 | sort -u > mrna_ids.txt
grep '^>' query_proteins.fasta | sed 's/^>//' | sort -u > query_ids.txt
n_query=$(wc -l < query_ids.txt)
n_match=$(comm -12 query_ids.txt mrna_ids.txt | wc -l)
if [[ "$n_query" -eq 0 || $((n_match * 100 / n_query)) -lt 95 ]]; then
  echo "Only ${n_match}/${n_query} protein headers match mRNA IDs of the GFF3" >&2
  exit 1
fi

# 2. reference proteome
cat $reference \
  | awk -v os="$exclude_os" '/^>/{keep = (os == "" || index($0, "OS=" os " ") == 0)} keep' \
  | sed '/^>/! s/[*.-]//g' > reference_proteins.fasta
grep -q '^>' reference_proteins.fasta || { echo "Empty reference proteome" >&2; exit 1; }

# 3. long reads -> one mRNA per primary alignment
isoseq_gff=""
if [[ -n "$bams" ]]; then
  : > longreads.gff3
  i=0
  for bam in $bams; do
    [[ -s "$bam" ]] || continue
    i=$((i + 1))
    samtools view -b -F 0x904 -q 20 "$bam" \
      | bedtools bamtobed -i stdin \
      | awk -v OFS='\t' -v s="lr${i}" '{print $1, "longread", "mRNA", $2 + 1, $3, ".", $6, ".", "ID=" s "_" NR ";"}' \
      >> longreads.gff3
  done
  [[ -s longreads.gff3 ]] && isoseq_gff="longreads.gff3"
fi

# 4. Mender config and run
cat > mender.cfg <<EOF
[input]
gff = coding_genes.gff3
genome_fa = ${genome}
proteome_fa = query_proteins.fasta
subject_fa = reference_proteins.fasta
isoseq_gff = ${isoseq_gff}
[output]
run_prefix = titan
results_dir = results
workdir = mender_workdir
db_dir = mender_db
[diamond]
threads = ${threads}
evalue = ${evalue}
[ids]
gene_template = ${gene_template}
trans_template = ${trans_template}
[merge_filters]
fix_partial = yes
skip_flags = ${skip_flags}
spanning_rescue = no
flags = ${flags}
[validation]
run_gt = yes
run_translation_validation = yes
final_gff_include_review = ${include_review}
run_agat = no
EOF
perl "$mender_dir/run_mender.pl" --config mender.cfg --steps 1,2,3,4,5,6,7 2>&1 | tee mender.log

# Mender builds an invalid model (overlapping exons / CDS) when the source genes overlap each other, and
# gffread then aborts in step 8: such merges are reverted to their source genes before validation.
if [[ -s results/titan/merges.gff ]]; then
  python3 - results/titan/merges.gff results/titan/removed.gff <<'PY'
import re, sys
merges_path, removed_path = sys.argv[1:3]

def gene_blocks(path):
    blocks = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) != 9:
            continue
        if f[2] == "gene":
            blocks.append([re.search(r"ID=([^;]+)", f[8]).group(1), f[8], [line]])
        elif blocks:
            blocks[-1][2].append(line)
    return blocks

def overlapping(lines):
    by_parent = {}
    for line in lines:
        f = line.split("\t")
        if f[2] in ("exon", "CDS"):
            parent = re.search(r"Parent=([^;\n]+)", f[8]).group(1)
            by_parent.setdefault((parent, f[2]), []).append((int(f[3]), int(f[4])))
    for iv in by_parent.values():
        iv.sort()
        if any(b[0] <= a[1] for a, b in zip(iv, iv[1:])):
            return True
    return False

merges = gene_blocks(merges_path)
removed = {b[0]: b for b in gene_blocks(removed_path)}
invalid = [b for b in merges if "merged_from=" in b[1] and overlapping(b[2])]
restore = set()
with open("mender_reverted_merges.tsv", "w") as out:
    out.write("merged_gene\tmerged_from\treason\n")
    for b in invalid:
        sources = re.search(r"merged_from=([^;]+)", b[1]).group(1)
        restore.update(sources.split(","))
        out.write(f"{b[0]}\t{sources}\toverlapping_exons_or_cds\n")
bad = {b[0] for b in invalid}
with open(merges_path, "w") as out:
    out.write("##gff-version 3\n")
    for b in [b for b in merges if b[0] not in bad] + [removed[g] for g in sorted(restore) if g in removed]:
        out.writelines(b[2])
with open(removed_path, "w") as out:
    out.write("##gff-version 3\n")
    for g, b in removed.items():
        if g not in restore:
            out.writelines(b[2])
print(f"Reverted {len(invalid)} merges with overlapping exons/CDS ({len(restore)} source genes restored)")
PY
fi
perl "$mender_dir/run_mender.pl" --config mender.cfg --steps 8,10 2>&1 | tee -a mender.log

validated="results/titan/validated.gff"
if [[ ! -s "$validated" ]]; then
  # no candidate passed the filters: Mender stops before writing a validated GFF
  echo "No validated GFF from Mender; annotation left unchanged" >&2
  cp "$gff3" final_annotation.gff3
  cp "$proteins_all" final_annotation_proteins_all.fasta
  cp "$proteins_main" final_annotation_proteins_main.fasta
  printf 'merged_gene\tchrom\tstart\tend\tstrand\tmerged_from\n' > mender_merges.tsv
  printf '{"merged_genes": 0, "source_genes_removed": 0}\n' > mender_summary.json
  exit 0
fi

# 5. TITAN-format outputs: input GFF3 minus the merged source genes, plus the merged genes, in genome order
python3 - "$validated" "$gff3" <<'PY'
import re, sys
validated, original = sys.argv[1:3]

def gene_blocks(path):
    blocks = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) != 9:
            continue
        if f[2] == "gene":
            blocks.append([f[0], int(f[3]), re.search(r"ID=([^;]+)", f[8]).group(1), f[8], [line]])
        elif blocks:
            blocks[-1][4].append(line)
    return blocks

original_blocks = gene_blocks(original)
chrom_order = {}
for b in original_blocks:
    chrom_order.setdefault(b[0], len(chrom_order))
merged = [b for b in gene_blocks(validated) if "merged_from=" in b[3]]

# Mender keeps the source CDS phases; recompute them along each merged transcript (first CDS kept)
def fix_phases(lines):
    cds = {}
    for i, line in enumerate(lines):
        f = line.rstrip("\n").split("\t")
        if f[2] == "CDS":
            cds.setdefault(re.search(r"Parent=([^;]+)", f[8]).group(1), []).append(i)
    for idx in cds.values():
        rows = [lines[i].rstrip("\n").split("\t") for i in idx]
        order = sorted(range(len(rows)), key=lambda k: int(rows[k][3]), reverse=rows[0][6] == "-")
        phase = int(rows[order[0]][7]) if rows[order[0]][7] in "012" else 0
        for k in order:
            rows[k][7] = str(phase)
            length = int(rows[k][4]) - int(rows[k][3]) + 1
            phase = (3 - (length - phase) % 3) % 3
        for k, i in enumerate(idx):
            lines[i] = "\t".join(rows[k]) + "\n"
    return lines

for b in merged:
    b[4] = fix_phases(b[4])
removed = {g for b in merged for g in re.search(r"merged_from=([^;]+)", b[3]).group(1).split(",")}
blocks = [b for b in original_blocks if b[2] not in removed] + merged
blocks.sort(key=lambda b: (chrom_order.get(b[0], len(chrom_order)), b[0], b[1]))
with open("final_annotation.gff3", "w") as out:
    out.write("##gff-version 3\n")
    for b in blocks:
        out.writelines(b[4] + ["###\n"])
PY
awk -F'\t' -v OFS='\t' 'BEGIN{print "merged_gene", "chrom", "start", "end", "strand", "merged_from"}
  $3=="gene" && $9 ~ /merged_from=/ {
    match($9, /ID=[^;]+/); id = substr($9, RSTART+3, RLENGTH-3)
    match($9, /merged_from=[^;]+/); mf = substr($9, RSTART+12, RLENGTH-12)
    print id, $1, $4, $5, $7, mf }' final_annotation.gff3 > mender_merges.tsv

awk -F'\t' 'NR>1{print $1}' mender_merges.tsv > merged_genes.txt
# merged genes with their children, by gene ID
python3 - <<'PY'
import re
merged = {l.strip() for l in open("merged_genes.txt") if l.strip()}
gene_of, out = {}, []
for line in open("final_annotation.gff3"):
    f = line.rstrip("\n").split("\t")
    if len(f) != 9:
        continue
    ident = re.search(r"(?:^|;)ID=([^;]+)", f[8])
    parent = re.search(r"(?:^|;)Parent=([^;,]+)", f[8])
    if f[2] == "gene":
        g = ident.group(1) if ident else None
    else:
        g = gene_of.get(parent.group(1)) if parent else None
        if g is None and parent and parent.group(1) in merged:
            g = parent.group(1)
    if ident and g:
        gene_of[ident.group(1)] = g
    if g in merged:
        out.append(line)
open("merged_genes.gff3", "w").writelines(out)
PY
if [[ -s merged_genes.gff3 ]]; then
  gffread -y merged_proteins.raw.fasta -g "$genome" merged_genes.gff3
else
  : > merged_proteins.raw.fasta
fi
python3 - "$proteins_all" "$proteins_main" <<'PY'
import re, sys, json
proteins_all, proteins_main = sys.argv[1:3]

def read_fasta(path):
    recs, name = [], None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs.append([name, []])
        elif name is not None:
            recs[-1][1].append(line)
    return [(n, "".join(s)) for n, s in recs]

transcript_gene = {}
for line in open("final_annotation.gff3"):
    f = line.rstrip("\n").split("\t")
    if len(f) == 9 and f[2] == "mRNA":
        t = re.search(r"ID=([^;]+)", f[8]).group(1)
        transcript_gene[t] = re.search(r"Parent=([^;,]+)", f[8]).group(1)

def tid(header):
    return re.sub(r"_CDS\d+\.prot$", "", header)

merged = {}
for name, seq in read_fasta("merged_proteins.raw.fasta"):
    seq = seq.replace(".", "*")
    merged.setdefault(transcript_gene.get(name), []).append((name + "_CDS1.prot", seq))
merged.pop(None, None)

def write(path_in, path_out, main):
    with open(path_out, "w") as out:
        for name, seq in read_fasta(path_in):
            if tid(name) in transcript_gene:
                out.write(f">{name}\n{seq}\n")
        for gene, recs in sorted(merged.items()):
            keep = [max(recs, key=lambda r: len(r[1]))] if main else recs
            for name, seq in keep:
                out.write(f">{name}\n{seq}\n")

write(proteins_all, "final_annotation_proteins_all.fasta", False)
write(proteins_main, "final_annotation_proteins_main.fasta", True)
removed = sum(len(l.rstrip("\n").split("\t")[5].split(",")) for l in list(open("mender_merges.tsv"))[1:])
json.dump({"merged_genes": len(merged), "source_genes_removed": removed}, open("mender_summary.json", "w"), indent=2)
PY
test -s final_annotation.gff3
test -s final_annotation_proteins_main.fasta
