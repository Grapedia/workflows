#!/usr/bin/env python3
"""Trace the reference genes transferred by Liftoff that are absent from the final annotation.

NOT part of the pipeline: an analysis tool. A Liftoff gene is "lost" when no CDS of the final annotation
overlaps its CDS on the same strand. Each lost gene is assigned to the first step where it disappears:

  prepare:<reason>       none of its transcripts left mikado prepare; the reason comes from prepare.log
                         ("internal_stops_N_mask" when the stops only exist on the hard-masked genome,
                         i.e. the soft-masked CDS has none, "internal_stops" otherwise)
  pick:in_host_utr       mikado pick kept no CDS over it and its CDS lies in the UTR of a coding Mikado
                         model (the host); host_gene, host_carried and host_origin describe that host
  pick:other             mikado pick kept no CDS over it for another reason
  finalize               a Mikado CDS covers it but no gene of the pre-filter annotation does
  filter:<reason>        the pre-filter gene covering it was removed by the TE / ncRNA filters
                         (genes_removed.tsv); that gene is not the Liftoff gene unless carried over

Outputs in --out-dir: lost_liftoff_genes.tsv (one row per lost gene) and lost_liftoff_summary.tsv
(per category: genes, expression and function in the reference, CDS hard-masked fraction).

Usage
  trace_lost_liftoff_genes.py --liftoff liftoff_previous_annotations.gff3 --final final_annotation.noTE.gff3 \
      --prefilter primary/final_annotation.gff3 --mikado-loci mikado.loci.gff3 \
      --mikado-prepared mikado_prepared.gtf --prepare-log prepare.log --genes-removed genes_removed.tsv \
      --correspondence liftoff_gene_id_correspondence.tsv --evidence-table gene_evidence_table.tsv \
      --hardmasked assembly_masked.EDTA.fasta --softmasked assembly_softmasked.EDTA.fasta \
      --ref-tpm gene_tpm_matrix.tsv --ref-func-dir functional/v43 --out-dir lost_liftoff
"""
import argparse
import collections
import os
import re
import subprocess
import sys
import tempfile

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_gene_evidence_table import MIN_TPM  # noqa: E402
from reference_qc_metrics import functional  # noqa: E402

COMP = str.maketrans("ACGTN", "TGCAN")


def attrs(field):
    return dict(kv.split("=", 1) for kv in field.strip().split(";") if "=" in kv)


def read_gff(path):
    """Genes, mRNA -> gene, CDS (chrom, start0, end, strand, phase) per mRNA, UTR intervals per mRNA."""
    genes, m2g = {}, {}
    cds, utr = collections.defaultdict(list), collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = attrs(f[8])
            if f[2] == "gene":
                genes[a["ID"]] = dict(chrom=f[0], start=int(f[3]), end=int(f[4]), strand=f[6], attrs=a)
            elif f[2] in ("mRNA", "transcript", "ncRNA", "lnc_RNA"):
                m2g[a["ID"]] = a.get("Parent")
            elif f[2] == "CDS":
                for p in a["Parent"].split(","):
                    cds[p].append((f[0], int(f[3]) - 1, int(f[4]), f[6], f[7]))
            elif "UTR" in f[2]:
                for p in a["Parent"].split(","):
                    utr[p].append((f[0], int(f[3]) - 1, int(f[4]), f[6]))
    return genes, m2g, cds, utr


def write_bed(path, rows):
    with open(path, "w") as out:
        for chrom, start, end, name, strand in sorted(rows):
            out.write("%s\t%d\t%d\t%s\t0\t%s\n" % (chrom, start, end, name, strand))


def cds_bed(path, m2g, cds):
    write_bed(path, [(c, s, e, m2g.get(t, t), st) for t, iv in cds.items() for c, s, e, st, _ in iv])


def overlaps(bedtools, a, b):
    """name in a -> set of names in b, same strand."""
    out = subprocess.run([bedtools, "intersect", "-a", a, "-b", b, "-s", "-wa", "-wb"], stdout=subprocess.PIPE,
                         stderr=subprocess.DEVNULL, universal_newlines=True, check=True).stdout
    hits = collections.defaultdict(set)
    for line in out.splitlines():
        x = line.split("\t")
        hits[x[3]].add(x[9])
    return hits


def load_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip().upper())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def cds_sequence(genome, intervals):
    iv = sorted(intervals, key=lambda x: x[1])
    seq = "".join(genome[c][s:e] for c, s, e, _, _ in iv)
    strand = iv[0][3]
    if strand == "-":
        seq = seq.translate(COMP)[::-1]
        phase = iv[-1][4]
    else:
        phase = iv[0][4]
    return seq[int(phase) if phase.isdigit() else 0:]


def internal_stops(seq):
    return sum(seq[i:i + 3] in ("TAA", "TAG", "TGA") for i in range(0, len(seq) - 5, 3))


def prepare_reasons(log_path):
    reasons = collections.defaultdict(set)
    pat = re.compile(r"Discarded generically invalid transcript liftoff_(\S+?)_t\d+, exception: (.*)")
    with open(log_path) as fh:
        for line in fh:
            m = pat.search(line)
            if not m:
                continue
            msg = m.group(2)
            reasons[m.group(1)].add("internal_stops" if "internal stop" in msg else
                                    "truncated_orf" if "truncated ORF" in msg else
                                    "utr_without_cds" if "no CDS" in msg else "other")
    return reasons


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for opt in ("liftoff", "final", "prefilter", "mikado-loci", "mikado-prepared", "prepare-log", "genes-removed",
                "correspondence", "evidence-table", "hardmasked", "softmasked", "ref-tpm", "ref-func-dir", "out-dir"):
        ap.add_argument("--" + opt, required=True)
    ap.add_argument("--bedtools", default="bedtools")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    tmp = tempfile.mkdtemp(dir=a.out_dir)

    lo_genes, lo_m2g, lo_cds, _ = read_gff(a.liftoff)
    lo_bed = os.path.join(tmp, "liftoff_cds.bed")
    cds_bed(lo_bed, lo_m2g, lo_cds)
    beds = {}
    for name, path in (("final", a.final), ("prefilter", a.prefilter), ("loci", a.mikado_loci)):
        _, m2g, cds, utr = read_gff(path)
        beds[name] = os.path.join(tmp, name + "_cds.bed")
        cds_bed(beds[name], m2g, cds)
        if name == "loci":
            beds["loci_utr"] = os.path.join(tmp, "loci_utr.bed")
            write_bed(beds["loci_utr"], [(c, s, e, m2g.get(t, t), st) for t, iv in utr.items() if t in cds
                                         for c, s, e, st in iv])
    hits = {k: overlaps(a.bedtools, lo_bed, beds[k]) for k in beds}

    coding = {lo_m2g[t] for t in lo_cds}
    lost = sorted(g for g in coding if g not in hits["final"])
    with open(a.mikado_prepared) as fh:
        text = fh.read()
    prepared = set(re.findall(r'transcript_id "liftoff_([^"]+)_t\d+"', text)) | set(re.findall(r'old_name "([^"]+)"', text))
    del text
    reasons = prepare_reasons(a.prepare_log)
    removed = pd.read_csv(a.genes_removed, sep="\t", index_col=0)
    corr = pd.read_csv(a.correspondence, sep="\t")
    carried = set(corr.loc[corr["decision"] == "carried_over_from_liftoff", "final_gene_id"])
    origin = pd.read_csv(a.evidence_table, sep="\t", index_col=0, usecols=["gene_id", "origin"])["origin"]

    hard, soft = load_fasta(a.hardmasked), load_fasta(a.softmasked)
    by_gene = collections.defaultdict(list)
    for t, iv in lo_cds.items():
        by_gene[lo_m2g[t]].append(iv)

    rows = []
    for g in lost:
        info = lo_genes[g]
        n_frac = max(cds_sequence(hard, iv).count("N") / max(sum(e - s for _, s, e, _, _ in iv), 1) for iv in by_gene[g])
        stops_soft = min(internal_stops(cds_sequence(soft, iv)) for iv in by_gene[g])
        row = dict(gene=g, chrom=info["chrom"], start=info["start"], end=info["end"], strand=info["strand"],
                   liftoff_coverage=float(info["attrs"].get("coverage", "nan")),
                   cds_hardmasked_frac=round(n_frac, 3), internal_stops_softmasked=stops_soft,
                   host_gene="", host_carried="", host_origin="", prefilter_genes="", removed_reason="")
        if g not in prepared:
            r = reasons.get(g, {"not_logged"})
            if "internal_stops" in r:
                r = {"internal_stops_N_mask" if stops_soft == 0 else "internal_stops"}
            row["step"] = "prepare:" + ",".join(sorted(r))
        elif g not in hits["loci"]:
            in_utr = g in hits["loci_utr"]
            row["step"] = "pick:in_host_utr" if in_utr else "pick:other"
            if in_utr:
                row["host_gene"] = ",".join(sorted(hits["loci_utr"][g]))
        elif g not in hits["prefilter"]:
            row["step"] = "finalize"
        else:
            pre = sorted(hits["prefilter"][g])
            row["prefilter_genes"] = ",".join(pre)
            rr = sorted({str(removed["reason"].get(x)) for x in pre if x in removed.index})
            row["removed_reason"] = ",".join(rr)
            row["step"] = "filter:" + (",".join(rr) or "unknown")
        rows.append(row)
    df = pd.DataFrame(rows).set_index("gene")

    # Hosts: map the Mikado host model to the pre-filter gene covering the same CDS (AEGIS renames genes)
    host_rows = [(g, h) for g in df.index for h in df.at[g, "host_gene"].split(",") if h]
    if host_rows:
        _, loci_m2g2, loci_cds, _ = read_gff(a.mikado_loci)
        hb = os.path.join(tmp, "host_cds.bed")
        hosts = {h for _, h in host_rows}
        write_bed(hb, [(c, s, e, loci_m2g2.get(t, t), st) for t, iv in loci_cds.items()
                       if loci_m2g2.get(t, t) in hosts for c, s, e, st, _ in iv])
        to_pre = overlaps(a.bedtools, hb, beds["prefilter"])
        for g in df.index[df["step"] == "pick:in_host_utr"]:
            pre = sorted({p for h in df.at[g, "host_gene"].split(",") for p in to_pre.get(h, ())})
            df.at[g, "host_carried"] = str(any(p in carried for p in pre))
            df.at[g, "host_origin"] = ",".join(sorted({str(origin.get(p)) for p in pre}))

    tpm = pd.read_csv(a.ref_tpm, sep="\t", index_col=0)
    df["ref_max_tpm"] = tpm.max(axis=1).reindex(df.index)
    df["ref_expressed"] = df["ref_max_tpm"] >= MIN_TPM
    df["ref_function"] = df.index.isin(functional(a.ref_func_dir)["any"])
    df.to_csv(os.path.join(a.out_dir, "lost_liftoff_genes.tsv"), sep="\t")

    summary = df.groupby("step").agg(
        genes=("chrom", "size"),
        ref_expressed_pct=("ref_expressed", lambda s: round(100 * s.mean(), 1)),
        ref_median_max_tpm=("ref_max_tpm", "median"),
        ref_function_pct=("ref_function", lambda s: round(100 * s.mean(), 1)),
        median_cds_hardmasked_frac=("cds_hardmasked_frac", "median"),
        liftoff_coverage_lt1=("liftoff_coverage", lambda s: int((s < 1).sum())))
    summary.to_csv(os.path.join(a.out_dir, "lost_liftoff_summary.tsv"), sep="\t")
    print("Liftoff coding genes: %d; lost: %d" % (len(coding), len(df)))
    print(summary.to_string())
    hosts = df[df["step"] == "pick:in_host_utr"]
    if len(hosts):
        print("in-host-UTR genes whose host is a carried-over gene: %d; host origins: %s" % (
            (hosts["host_carried"] == "True").sum(),
            hosts["host_origin"].str.split(",").explode().value_counts().to_dict()))
    for f in os.listdir(tmp):
        os.remove(os.path.join(tmp, f))
    os.rmdir(tmp)


if __name__ == "__main__":
    main()
