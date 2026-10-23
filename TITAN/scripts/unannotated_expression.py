#!/usr/bin/env python3
"""Expressed loci that an annotation misses ("expression outside the annotation").

NOT part of the pipeline: an analysis tool. Needs bedtools on PATH.

Loci = genes of the multi-sample StringTie merges (first merge, then loci of the next merges that do not overlap
an earlier one). A locus is expressed when transcripts assembled in >= --min-samples samples with
TPM >= --min-tpm (per-sample StringTie GTFs) overlap its exons, and its merged exons span >= --min-exon-bp.
A locus is missed by an annotation when none of its exonic bases overlaps an exon of that annotation (strand
ignored: part of the RNA-seq is unstranded). Loci missed by the first annotation are then labelled, in this
order, with the first matching explanation track (--explain name=gff, exon/feature overlap; the special
--explain te uses >= --te-frac of the exonic bases in --te); the rest are "unexplained".

Outputs in --out-dir: loci.tsv (one row per expressed locus) and summary.json.
"""
import argparse
import collections
import json
import os
import re
import subprocess
import tempfile


def run(cmd, out):
    with open(out, "w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)


def gtf_attr(attrs, key):
    m = re.search(key + r' "([^"]+)"', attrs)
    return m.group(1) if m else None


def exon_bed_from_gff(path, out, types=("exon",)):
    """Exon-like features of a GFF3/GTF as BED (falls back to CDS/gene when a feature has no exons)."""
    with open(path) as fh, open(out, "w") as o:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == 9 and f[2] in types:
                o.write(f"{f[0]}\t{int(f[3]) - 1}\t{f[4]}\n")


def sort_bed(path):
    run(["bedtools", "sort", "-i", path], path + ".sorted")
    os.replace(path + ".sorted", path)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged-gtfs", nargs="+", required=True, help="StringTie --merge GTFs, by priority")
    ap.add_argument("--sample-gtfs", required=True, help="file listing the per-sample StringTie GTFs")
    ap.add_argument("--annotation", action="append", required=True, help="name=gff3; the first one is explained")
    ap.add_argument("--explain", action="append", default=[], help="name=gff, in priority order")
    ap.add_argument("--explain-types", default="exon,CDS,tRNA,ncRNA,lncRNA,transcript,repeat_region",
                    help="feature types used from explanation tracks")
    ap.add_argument("--te", help="TE annotation (EDTA TEanno.gff3)")
    ap.add_argument("--te-frac", type=float, default=0.5)
    ap.add_argument("--min-tpm", type=float, default=1.0)
    ap.add_argument("--min-samples", type=int, default=3)
    ap.add_argument("--min-exon-bp", type=int, default=200)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    tmp = tempfile.mkdtemp(dir=args.out_dir, prefix="tmp_")

    # 1. loci from the merges (exons as BED, name = locus)
    locus_exons = os.path.join(tmp, "locus_exons.bed")
    taken = os.path.join(tmp, "taken.bed")
    open(taken, "w").close()
    n_exons = collections.Counter()
    with open(locus_exons, "w") as out:
        for k, gtf in enumerate(args.merged_gtfs):
            rows = []
            for line in open(gtf):
                f = line.rstrip("\n").split("\t")
                if len(f) == 9 and f[2] == "exon":
                    tid = gtf_attr(f[8], "transcript_id")
                    rows.append((f[0], int(f[3]) - 1, int(f[4]), f"m{k}:{gtf_attr(f[8], 'gene_id')}", tid))
            cand = os.path.join(tmp, f"cand{k}.bed")
            with open(cand, "w") as o:
                for r in rows:
                    o.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")
            sort_bed(cand)
            overl = set()
            if os.path.getsize(taken):
                res = subprocess.run(["bedtools", "intersect", "-a", cand, "-b", taken, "-u"],
                                     capture_output=True, text=True, check=True).stdout
                overl = {l.split("\t")[3] for l in res.splitlines()}
            per_tx = collections.Counter(r[4] for r in rows)
            for r in rows:
                if r[3] not in overl:
                    out.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")
                    n_exons[r[3]] = max(n_exons[r[3]], per_tx[r[4]])
            with open(taken, "a") as t:
                for r in rows:
                    if r[3] not in overl:
                        t.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")
            sort_bed(taken)
    sort_bed(locus_exons)
    # merged exons and exonic bp per locus
    loc_bp, loc_pos = collections.Counter(), {}
    lm = os.path.join(tmp, "locus_merged.bed")
    per_locus = collections.defaultdict(list)
    for line in open(locus_exons):
        c, s, e, name = line.rstrip("\n").split("\t")
        per_locus[name].append((c, int(s), int(e)))
    with open(lm, "w") as o:
        for name, iv in per_locus.items():
            iv.sort()
            cur = None
            for c, s, e in iv:
                if cur and s <= cur[2]:
                    cur[2] = max(cur[2], e)
                else:
                    if cur:
                        o.write(f"{cur[0]}\t{cur[1]}\t{cur[2]}\t{name}\n"); loc_bp[name] += cur[2] - cur[1]
                    cur = [c, s, e]
            o.write(f"{cur[0]}\t{cur[1]}\t{cur[2]}\t{name}\n"); loc_bp[name] += cur[2] - cur[1]
            loc_pos[name] = (iv[0][0], min(x[1] for x in iv), max(x[2] for x in iv))
    sort_bed(lm)

    # 2. sample support
    samp = os.path.join(tmp, "sample_exons.bed")
    with open(samp, "w") as o:
        for gtf in [l.strip() for l in open(args.sample_gtfs) if l.strip()]:
            sname = os.path.basename(gtf).split("_transcriptome")[0]
            tpm = {}
            for line in open(gtf):
                f = line.rstrip("\n").split("\t")
                if len(f) == 9 and f[2] == "transcript":
                    t = gtf_attr(f[8], "TPM")
                    tpm[gtf_attr(f[8], "transcript_id")] = float(t) if t else 0.0
            for line in open(gtf):
                f = line.rstrip("\n").split("\t")
                if len(f) == 9 and f[2] == "exon":
                    v = tpm.get(gtf_attr(f[8], "transcript_id"), 0.0)
                    if v >= args.min_tpm:
                        o.write(f"{f[0]}\t{int(f[3]) - 1}\t{f[4]}\t{sname}\t{v}\n")
    sort_bed(samp)
    res = subprocess.run(["bedtools", "intersect", "-a", lm, "-b", samp, "-wa", "-wb", "-sorted"],
                         capture_output=True, text=True, check=True).stdout
    support, maxtpm = collections.defaultdict(set), collections.Counter()
    for line in res.splitlines():
        f = line.split("\t")
        support[f[3]].add(f[7])
        maxtpm[f[3]] = max(maxtpm[f[3]], float(f[8]))
    expressed = [n for n in per_locus if len(support[n]) >= args.min_samples and loc_bp[n] >= args.min_exon_bp]
    ex_bed = os.path.join(tmp, "expressed.bed")
    with open(ex_bed, "w") as o:
        keep = set(expressed)
        for line in open(lm):
            if line.rstrip("\n").split("\t")[3] in keep:
                o.write(line)

    # 3. overlap with annotations and explanation tracks
    def covered(track_bed):
        out = subprocess.run(["bedtools", "intersect", "-a", ex_bed, "-b", track_bed, "-u"],
                             capture_output=True, text=True, check=True).stdout
        return {l.split("\t")[3] for l in out.splitlines()}

    ann_names, hits = [], {}
    for spec in args.annotation:
        name, path = spec.split("=", 1)
        bed = os.path.join(tmp, f"ann_{name}.bed")
        exon_bed_from_gff(path, bed)
        sort_bed(bed)
        hits[name] = covered(bed)
        ann_names.append(name)
    types = tuple(args.explain_types.split(","))
    te_frac = {}
    if args.te:
        te_bed = os.path.join(tmp, "te.bed")
        with open(args.te) as fh, open(te_bed, "w") as o:
            for line in fh:
                f = line.split("\t")
                if len(f) == 9 and not line.startswith("#"):
                    o.write(f"{f[0]}\t{int(f[3]) - 1}\t{f[4]}\n")
        sort_bed(te_bed)
        run(["bedtools", "merge", "-i", te_bed], te_bed + ".m")
        out = subprocess.run(["bedtools", "intersect", "-a", ex_bed, "-b", te_bed + ".m", "-wo"],
                             capture_output=True, text=True, check=True).stdout
        bp = collections.Counter()
        for line in out.splitlines():
            f = line.split("\t")
            bp[f[3]] += int(f[-1])
        te_frac = {n: bp[n] / loc_bp[n] for n in expressed}
    explain = []
    for spec in args.explain:
        if spec == "te":
            explain.append(("te", {n for n in expressed if te_frac.get(n, 0) >= args.te_frac}))
            continue
        name, path = spec.split("=", 1)
        bed = os.path.join(tmp, f"exp_{name}.bed")
        exon_bed_from_gff(path, bed, types)
        sort_bed(bed)
        explain.append((name, covered(bed)))

    first = ann_names[0]
    rows, why = [], collections.Counter()
    for n in expressed:
        missed = {a: n not in hits[a] for a in ann_names}
        reason = ""
        if missed[first]:
            for name, hset in explain:
                if n in hset:
                    reason = name
                    break
            reason = reason or ("unexplained_multiexon" if n_exons[n] >= 2 else "unexplained_monoexon")
            why[reason] += 1
        c, s, e = loc_pos[n]
        rows.append([n, c, s + 1, e, loc_bp[n], n_exons[n], len(support[n]), round(maxtpm[n], 2),
                     round(te_frac.get(n, 0), 2)] + [int(missed[a]) for a in ann_names] + [reason])
    with open(os.path.join(args.out_dir, "loci.tsv"), "w") as o:
        o.write("\t".join(["locus", "chrom", "start", "end", "exonic_bp", "max_exons", "samples", "max_tpm", "te_frac"]
                          + [f"missed_by_{a}" for a in ann_names] + [f"reason_missed_by_{first}"]) + "\n")
        for r in rows:
            o.write("\t".join(map(str, r)) + "\n")
    summary = {"params": vars(args), "loci_total": len(per_locus), "loci_expressed": len(expressed),
               "missed": {a: sum(1 for n in expressed if n not in hits[a]) for a in ann_names},
               f"reason_missed_by_{first}": dict(why)}
    for a in ann_names[1:]:
        summary[f"missed_by_{first}_and_{a}"] = sum(1 for n in expressed if n not in hits[first] and n not in hits[a])
        summary[f"missed_by_{first}_only_vs_{a}"] = sum(1 for n in expressed if n not in hits[first] and n in hits[a])
        summary[f"missed_by_{a}_only_vs_{first}"] = sum(1 for n in expressed if n in hits[first] and n not in hits[a])
    json.dump(summary, open(os.path.join(args.out_dir, "summary.json"), "w"), indent=2)
    subprocess.run(["rm", "-rf", tmp])
    print(json.dumps({k: v for k, v in summary.items() if k != "params"}, indent=2))


if __name__ == "__main__":
    main()
