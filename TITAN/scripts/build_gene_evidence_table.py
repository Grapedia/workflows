#!/usr/bin/env python3
"""Build a one-row-per-gene evidence + structure table for a finished TITAN run.

Purpose: characterise genes of final_annotation.gff3 that carry no independent
support (no RNA-seq expression, no functional hit) so they can be studied and,
if a rule can be justified, set aside.  See docs/user/unsupported_genes.md.

Inputs are all read from a TITAN production output directory (the 5 numbered
folders).  Output: a TSV, one row per gene.  Python >= 3.6, pandas required.

Usage:
    build_gene_evidence_table.py --outdir data/titan_prod_out \
        --out data/unsupported_genes_analysis/gene_evidence_table.tsv
"""
import argparse
import collections
import math
import os
import re
import subprocess
import sys
import tempfile

import pandas as pd

MIN_TPM = 0.5
# InterProScan analyses that describe sequence composition, not function/family.
NON_INFORMATIVE_IPR = {"MobiDBLite", "Coils"}
# Text of functional hits that identifies a transposable-element-encoded protein.
TE_TEXT = re.compile(
    r"transposon|transposase|retrotranspos|reverse transcriptase|\bgag\b|gag-pol|polyprotein|"
    r"copia|gypsy|\bTy1\b|\bTy3\b|retrovir|retroelement|mutator|\bMULE\b|helitron|"
    r"hAT family|HAT, C-terminal|CACTA|\btc1\b|mariner|harbinger|piggybac|\bLINE-1\b|non-LTR|"
    r"integrase|\bDDE\b", re.I)
CDS_SUFFIX = re.compile(r"_t\d+_CDS\d+\.prot$")


def parse_attrs(field):
    out = {}
    for kv in field.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def read_gff(path):
    """Return genes (dict id -> info) with per-mRNA exon/CDS structure."""
    genes, mrna_parent, mrnas = {}, {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, src, typ, s, e, score, strand = f[0], f[1], f[2], int(f[3]), int(f[4]), f[5], f[6]
            a = parse_attrs(f[8])
            if typ == "gene":
                genes[a["ID"]] = dict(chrom=chrom, start=s, end=e, strand=strand)
            elif typ in ("mRNA", "ncRNA"):
                mrna_parent[a["ID"]] = a["Parent"]
                mrnas[a["ID"]] = dict(gene=a["Parent"], type=typ, start=s, end=e,
                                      exons=[], cds=[])
            elif typ == "exon":
                for p in a["Parent"].split(","):
                    mrnas[p]["exons"].append((s, e))
            elif typ == "CDS":
                for p in a["Parent"].split(","):
                    mrnas[p]["cds"].append((s, e))
    per_gene = collections.defaultdict(list)
    for mid, m in mrnas.items():
        per_gene[m["gene"]].append((mid, m))
    rows = {}
    for gid, g in genes.items():
        ms = per_gene.get(gid, [])
        if not ms:
            continue
        # main isoform = longest CDS (ties -> most exons)
        def key(item):
            m = item[1]
            return (sum(e - s + 1 for s, e in m["cds"]), len(m["exons"]))
        mid, m = max(ms, key=key)
        cds_len = sum(e - s + 1 for s, e in m["cds"])
        exon_len = sum(e - s + 1 for s, e in m["exons"])
        rows[gid] = dict(
            g,
            main_mrna=mid,
            biotype=m["type"],
            n_mrna=len(ms),
            n_exons=len(m["exons"]),
            n_cds_exons=len(m["cds"]),
            gene_len=g["end"] - g["start"] + 1,
            exon_len=exon_len,
            cds_len=cds_len,
            utr_len=max(exon_len - cds_len, 0),
            single_exon=len(m["exons"]) == 1,
            cds_intervals=sorted(m["cds"]),
        )
    return rows


def read_mikado(path):
    """Primary-mRNA attributes of Mikado loci, keyed by (chrom, start, end, strand)."""
    genes, out = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            typ = f[2]
            a = parse_attrs(f[8])
            if typ == "gene":
                genes[a["ID"]] = (f[0], int(f[3]), int(f[4]), f[6])
            elif typ == "mRNA":
                gid = a["Parent"]
                if gid not in genes:
                    continue
                if a.get("primary", "False") != "True" and genes[gid] in out:
                    continue
                alias = a.get("alias", "")
                out[genes[gid]] = dict(
                    mikado_alias=alias,
                    mikado_evidence=a.get("evidence"),
                    mikado_count=a.get("count"),
                    mikado_canon_prop=a.get("canonical_proportion"),
                    mikado_has_start=a.get("has_start_codon"),
                    mikado_has_stop=a.get("has_stop_codon"),
                    mikado_is_reference=a.get("is_reference"),
                    mikado_retained_intron=a.get("retained_intron"),
                    mikado_ccode=a.get("ccode"),
                    mikado_multiexonic=None,
                )
    return out


def origin_from_alias(alias, is_reference):
    """Collapse Mikado's `alias` (winning-model name) to the evidence family."""
    a = alias.lower()
    if a.startswith("braker") or a.startswith("augustus") or a.startswith("genemark"):
        return "ab_initio_braker"
    if a.startswith("helixer"):
        return "ab_initio_helixer"
    if a.startswith("egapx") or a.startswith("gnomon"):
        return "egapx_gnomon"
    if a.startswith("liftoff") or a.startswith("vit") or a.startswith("ugt") \
            or a.startswith("vitvi") or is_reference == "True":
        return "liftoff_previous"
    if a.startswith("star") or a.startswith("psiclass") or a.startswith("stringtie") \
            or a.startswith("long") or a.startswith("minimap"):
        return "rnaseq_assembly"
    return "other"


def gene_of_query(q):
    return CDS_SUFFIX.sub("", q)


def load_expression(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    return pd.DataFrame({
        "max_tpm": df.max(axis=1),
        "n_samples_expr": (df >= MIN_TPM).sum(axis=1),
        "n_samples_any": (df > 0).sum(axis=1),
        "mean_tpm": df.mean(axis=1),
    }), df.shape[1]


def load_functional(outdir):
    fa = os.path.join(outdir, "02_functional_annotation")
    fn = {}
    # Diamond2GO
    d = pd.read_csv(os.path.join(fa, "diamond2go", "final_annotation_proteins_main.diamond2go.tsv"),
                    sep="\t", dtype=str)
    d.columns = [c.lstrip("#") for c in d.columns]
    d["gene"] = d["gene_id"].map(gene_of_query)
    fn["d2go_any"] = set(d["gene"])
    fn["d2go_go"] = set(d.loc[d["GO-term"].notna() & (d["GO-term"] != "-"), "gene"])
    # eggNOG
    e = pd.read_csv(os.path.join(fa, "eggnog", "final_annotation_proteins_main.emapper.annotations"),
                    sep="\t", comment="#", header=None, dtype=str, low_memory=False)
    e.columns = ["query", "seed", "evalue", "score", "OGs", "max_lvl", "COG", "Description",
                 "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
                 "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG", "PFAMs"][:e.shape[1]]
    e["gene"] = e["query"].map(gene_of_query)
    fn["egg_any"] = set(e["gene"])
    fn["egg_desc"] = set(e.loc[e["Description"].notna() & ~e["Description"].isin(["-", ""]), "gene"])
    fn["egg_go"] = set(e.loc[e["GOs"].notna() & (e["GOs"] != "-"), "gene"])
    fn["egg_kegg"] = set(e.loc[e["KEGG_ko"].notna() & (e["KEGG_ko"] != "-"), "gene"])
    # InterProScan
    cols = ["query", "md5", "len", "analysis", "sig_acc", "sig_desc", "start", "end", "score",
            "status", "date", "ipr", "ipr_desc", "go", "pathways"]
    ip = pd.read_csv(os.path.join(fa, "interproscan", "final_annotation_proteins_main.tsv"),
                     sep="\t", header=None, names=cols, usecols=range(13), dtype=str,
                     low_memory=False)
    ip["gene"] = ip["query"].map(gene_of_query)
    fn["ipr_any"] = set(ip["gene"])
    inf = ip[~ip["analysis"].isin(NON_INFORMATIVE_IPR)]
    fn["ipr_informative"] = set(inf["gene"])
    fn["ipr_entry"] = set(ip.loc[ip["ipr"].notna() & (ip["ipr"] != "-"), "gene"])
    # TE-encoded protein according to any of the three functional methods
    ipr_txt = ip["sig_desc"].fillna("") + " | " + ip["ipr_desc"].fillna("")
    te = set(ip.loc[ipr_txt.map(lambda t: bool(TE_TEXT.search(t))), "gene"])
    te |= set(e.loc[e["Description"].fillna("").map(lambda t: bool(TE_TEXT.search(t))), "gene"])
    te |= set(d.loc[d["gene_name"].fillna("").map(lambda t: bool(TE_TEXT.search(t))), "gene"])
    fn["te_domain_hit"] = te
    # MapMan / Mercator4: IDs are quoted and lower-case; bin 99 = "not assigned",
    # bin 50 = "uncharacterised context" -> only bins 1-30 count as a functional call.
    mm = pd.read_csv(os.path.join(fa, "mercator4", "PNT2T.results.txt"), sep="\t", dtype=str)
    mm = mm[mm["TYPE"] == "T"].copy()
    mm["gene"] = mm["IDENTIFIER"].str.strip("'").str.lower()
    mm["top"] = mm["BINCODE"].str.strip("'").str.split(".").str[0].astype(int)
    fn["mapman_any"] = set(mm.loc[mm["top"] != 99, "gene"])
    fn["mapman_informative"] = set(mm.loc[mm["top"] <= 30, "gene"])
    return fn


def load_omamer(path):
    df = pd.read_csv(path, sep="\t", comment="!", dtype={"qseqid": str})
    df["gene"] = df["qseqid"].map(gene_of_query)
    df = df.sort_values("family_p", ascending=False).drop_duplicates("gene")
    return df.set_index("gene")[["hoglevel", "family_p", "family_count", "subfamily_score",
                                 "qseq_overlap"]].rename(columns=lambda c: "omamer_" + c)


def protein_features(path):
    """Length, low-complexity and composition of the main protein of each gene."""
    feats, name, seq = {}, None, []

    def flush():
        if name is None:
            return
        s = "".join(seq).replace("*", "").replace("-", "")
        if not s:
            return
        cnt = collections.Counter(s)
        n = len(s)
        ent = -sum(c / n * math.log(c / n, 2) for c in cnt.values())
        feats[gene_of_query(name)] = dict(
            prot_len=n,
            prot_entropy=round(ent, 3),
            prot_top_aa_frac=round(max(cnt.values()) / n, 3),
            prot_start_M=s[0] == "M",
            prot_frac_X=round(cnt.get("X", 0) / n, 4),
        )

    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    flush()
    return pd.DataFrame.from_dict(feats, orient="index")


def te_overlap(genes, gene_df, te_gff, tmp):
    """Fraction of CDS bp and of gene-span bp overlapped by an EDTA TE annotation."""
    te_bed = os.path.join(tmp, "te.bed")
    with open(te_gff) as fh, open(te_bed, "w") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) > 8 and f[2] not in ("repeat_region", "target_site_duplication"):
                out.write("%s\t%d\t%d\n" % (f[0], int(f[3]) - 1, int(f[4])))
    cds_bed, gene_bed = os.path.join(tmp, "cds.bed"), os.path.join(tmp, "gene.bed")
    with open(cds_bed, "w") as c, open(gene_bed, "w") as g:
        for gid, r in genes.items():
            g.write("%s\t%d\t%d\t%s\n" % (r["chrom"], r["start"] - 1, r["end"], gid))
            for s, e in r["cds_intervals"]:
                c.write("%s\t%d\t%d\t%s\n" % (r["chrom"], s - 1, e, gid))
    res = {}
    for label, bed in (("te_frac_cds", cds_bed), ("te_frac_gene", gene_bed)):
        sortd = os.path.join(tmp, label + ".sorted.bed")
        subprocess.check_call("sort -k1,1 -k2,2n %s > %s" % (bed, sortd), shell=True)
        te_sorted = os.path.join(tmp, "te.sorted.bed")
        subprocess.check_call("sort -k1,1 -k2,2n %s | bedtools merge -i - > %s" % (te_bed, te_sorted),
                              shell=True)
        out = subprocess.check_output(
            "bedtools intersect -a %s -b %s | awk -F'\\t' '{d[$4]+=$3-$2} END{for(k in d)print k\"\\t\"d[k]}'"
            % (sortd, te_sorted), shell=True).decode()
        ov = {l.split("\t")[0]: int(l.split("\t")[1]) for l in out.strip().split("\n") if l}
        denom = (gene_df["cds_len"] if label == "te_frac_cds" else gene_df["gene_len"])
        res[label] = pd.Series(ov).reindex(gene_df.index).fillna(0) / denom.clip(lower=1)
    return pd.DataFrame(res)


def neighbourhood(gene_df, genes):
    """Overlap with other genes (sense/antisense, CDS-level) and distance to neighbours."""
    rows = {}
    for chrom, sub in gene_df.groupby("chrom"):
        sub = sub.sort_values("start")
        ids = list(sub.index)
        starts, ends, strands = sub["start"].values, sub["end"].values, sub["strand"].values
        n = len(ids)
        max_end_so_far = -1
        for i, gid in enumerate(ids):
            left = starts[i] - ends[i - 1] - 1 if i > 0 else None
            right = starts[i + 1] - ends[i] - 1 if i + 1 < n else None
            rows[gid] = dict(dist_prev=left, dist_next=right)
        # overlap partners (sweep; genes are ~56k so O(n*k) is fine)
        j = 0
        for i in range(n):
            k = i + 1
            while k < n and starts[k] <= ends[i]:
                for a, b in ((i, k), (k, i)):
                    key = "ovl_same_strand" if strands[a] == strands[b] else "ovl_antisense"
                    rows[ids[a]][key] = rows[ids[a]].get(key, 0) + 1
                    # CDS-level competition, same strand only
                    if strands[a] == strands[b]:
                        ca = genes[ids[a]]["cds_intervals"]
                        cb = genes[ids[b]]["cds_intervals"]
                        if any(s1 <= e2 and s2 <= e1 for s1, e1 in ca for s2, e2 in cb):
                            rows[ids[a]]["ovl_cds_same_strand"] = rows[ids[a]].get("ovl_cds_same_strand", 0) + 1
                    # is a fully nested inside b ?
                    if starts[b] <= starts[a] and ends[a] <= ends[b]:
                        rows[ids[a]]["nested_in_other"] = 1
                k += 1
    df = pd.DataFrame.from_dict(rows, orient="index")
    for c in ("ovl_same_strand", "ovl_antisense", "ovl_cds_same_strand", "nested_in_other"):
        if c not in df:
            df[c] = 0
    return df.fillna({"ovl_same_strand": 0, "ovl_antisense": 0, "ovl_cds_same_strand": 0,
                      "nested_in_other": 0})


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True, help="TITAN output dir (holds 01_final_annotation ...)")
    ap.add_argument("--out", required=True, help="output TSV")
    a = ap.parse_args()
    o = a.outdir
    p = lambda *x: os.path.join(o, *x)

    print("parsing final GFF3 ...", file=sys.stderr)
    genes = read_gff(p("01_final_annotation", "primary", "final_annotation.gff3"))
    gene_df = pd.DataFrame.from_dict(genes, orient="index")
    gene_df.index.name = "gene_id"

    print("joining Mikado attributes (by coordinates) ...", file=sys.stderr)
    mk = read_mikado(p("04_evidence", "mikado", "final_mikado_annotation.gff3"))
    mk_rows = {gid: mk.get((r["chrom"], r["start"], r["end"], r["strand"]), {}) for gid, r in genes.items()}
    mk_df = pd.DataFrame.from_dict(mk_rows, orient="index")
    mk_df["mikado_matched"] = mk_df["mikado_alias"].notna()
    mk_df["origin"] = [origin_from_alias(x, y) if isinstance(x, str) else "unmatched"
                       for x, y in zip(mk_df["mikado_alias"], mk_df["mikado_is_reference"])]
    gene_df = gene_df.join(mk_df)

    print("expression ...", file=sys.stderr)
    expr, n_samples = load_expression(p("01_final_annotation", "quality_report", "expression_validation",
                                        "gene_tpm_matrix.tsv"))
    gene_df = gene_df.join(expr)
    gene_df["expressed"] = gene_df["max_tpm"] >= MIN_TPM

    print("functional annotation ...", file=sys.stderr)
    fn = load_functional(o)
    for k, v in fn.items():
        gene_df[k] = (gene_df.index.str.lower() if k.startswith("mapman") else gene_df.index).isin(v)
    gene_df["functional_any"] = gene_df[["d2go_any", "egg_any", "ipr_any"]].any(axis=1)
    gene_df["has_go"] = gene_df[["d2go_go", "egg_go"]].any(axis=1)
    gene_df["functional_strict"] = gene_df[["d2go_go", "egg_desc", "egg_go", "egg_kegg", "ipr_informative"]].any(axis=1)

    print("liftoff ...", file=sys.stderr)
    lo = pd.read_csv(p("01_final_annotation", "primary", "liftoff_gene_id_correspondence.tsv"), sep="\t")
    gene_df["liftoff_carried"] = gene_df.index.isin(lo.loc[lo["decision"] == "carried_over_from_liftoff", "final_gene_id"])

    print("OMAMER / proteins ...", file=sys.stderr)
    gene_df = gene_df.join(load_omamer(p("01_final_annotation", "quality_report", "omark", "proteins_main.omamer")))
    gene_df["omamer_hit"] = gene_df["omamer_hoglevel"].notna()
    gene_df = gene_df.join(protein_features(p("01_final_annotation", "primary", "final_annotation_proteins_main.fasta")))

    print("TE overlap ...", file=sys.stderr)
    with tempfile.TemporaryDirectory() as tmp:
        gene_df = gene_df.join(te_overlap(genes, gene_df, p("05_run_info", "intermediate_files",
                                                            "evidence_data", "EDTA", "edta.TEanno.gff3"), tmp))

    print("neighbourhood ...", file=sys.stderr)
    gene_df = gene_df.join(neighbourhood(gene_df, genes))

    gene_df["unsupported"] = ~gene_df["expressed"] & ~gene_df["functional_any"]
    gene_df["zero_evidence"] = gene_df["unsupported"] & ~gene_df["liftoff_carried"]
    gene_df = gene_df.drop(columns=["cds_intervals"])
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    gene_df.to_csv(a.out, sep="\t")
    print("wrote %s (%d genes, %d expression samples)" % (a.out, len(gene_df), n_samples), file=sys.stderr)


if __name__ == "__main__":
    main()
