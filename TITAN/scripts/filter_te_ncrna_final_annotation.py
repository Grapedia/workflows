#!/usr/bin/env python3
"""Build the protein-coding, TE-cleaned variant of the final annotation (primary_noTE).

Starting from final_annotation.gff3, this removes
  1. every non-coding gene (no mRNA child, i.e. ncRNA genes), and
  2. every TE-derived gene: >= --te-frac (default 0.80) of the CDS of its main isoform lies in an
     EDTA TE annotation AND the gene has no independent reason to be kept:
       - not carried over from the Liftoff reference annotation (liftoff_gene_id_correspondence.tsv),
       - >= 50 % of its CDS does not overlap (same strand) the Liftoff-transferred previous annotation,
       - not a BUSCO match,
       - no non-Vitis homolog (DIAMOND best hit >= 35 % identity and >= 50 % query coverage),
       - no host (non-TE) functional hit: Diamond2GO GO term, eggNOG description/GO/KEGG KO, or an
         informative InterProScan analysis (MobiDBLite/Coils do not count), unless that very
         protein also carries a TE-encoded domain/description.
  3. with --te-function-terms, every gene whose main protein carries one of the listed TE functional
     annotations (InterPro / Pfam / PANTHER accession, exact eggNOG or Diamond2GO description), whatever
     its other annotations, unless it is carried over from the Liftoff reference annotation.
It then filters the protein FASTAs and the Diamond2GO / eggNOG / InterProScan tables the same way.

The original files are never modified.  MapMan/Mercator4 is a manual step outside the pipeline and is
not handled here.  Stdlib only (runs in the plain python container).

A third sub-command, `compare`, tabulates the quality of several candidate annotations side by side
(gene counts, TE content, function, expression, BUSCO, OMArk): see candidates_comparison.tsv.

Two steps, because the homology search sits in between (it is run only on the few thousand candidates):
  select    -> candidates.tsv + candidates.fasta (TE-derived genes not protected by any other evidence)
  finalize  -> reads the DIAMOND hits of candidates.fasta, drops the candidates without a non-Vitis
               homolog, the genes with a listed TE function (step 3) and the ncRNA genes, and writes the
               filtered GFF3 / FASTAs / functional tables.
This is the in-pipeline version of scripts/filter_te_genes.py (post-hoc analysis tool).
"""
import argparse
import bisect
import csv
import json
import os
import re
from collections import defaultdict

# InterProScan analyses that describe sequence composition, not function/family.
NON_INFORMATIVE_IPR = {"MobiDBLite", "Coils"}
# Text of functional hits that identifies a transposable-element-encoded protein.
TE_TEXT = re.compile(
    r"transposon|transposase|retrotranspos|reverse transcriptase|\bgag\b|gag-pol|polyprotein|"
    r"copia|gypsy|\bTy1\b|\bTy3\b|retrovir|retroelement|mutator|\bMULE\b|helitron|"
    r"hAT family|HAT, C-terminal|CACTA|\btc1\b|mariner|harbinger|piggybac|\bLINE-1\b|non-LTR|"
    r"integrase|\bDDE\b", re.I)
# AEGIS protein headers: <gene>_t001_CDS1.prot
CDS_RECORD_SUFFIX = re.compile(r"_t\d+_CDS\d+(\.prot)?$")
# TE annotation features that are not the TE body itself
TE_SKIP_TYPES = {"repeat_region", "target_site_duplication"}


def parse_attrs(field):
    out = {}
    for item in field.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def gene_of(query):
    return CDS_RECORD_SUFFIX.sub("", query)


def read_gff3(path):
    """genes {id: (chrom,start,end)}, mrna -> gene, per-mRNA CDS intervals."""
    genes, mrna_gene, cds = {}, {}, defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = parse_attrs(f[8])
            if f[2] == "gene":
                genes[a["ID"]] = f[0]
            elif f[2] in ("mRNA", "transcript") and "Parent" in a:
                mrna_gene[a["ID"]] = a["Parent"].split(",")[0]
            elif f[2] == "CDS" and "Parent" in a:
                for p in a["Parent"].split(","):
                    cds[p].append((f[0], int(f[3]), int(f[4]), f[6]))
    return genes, mrna_gene, cds


def read_te_spans(path):
    raw = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 8 and f[2] not in TE_SKIP_TYPES:
                raw[f[0]].append((int(f[3]), int(f[4])))
    spans, starts = {}, {}
    for chrom, lst in raw.items():
        lst.sort()
        merged = []
        for s, e in lst:
            if merged and s <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        spans[chrom] = merged
        starts[chrom] = [s for s, _ in merged]
    return spans, starts


def overlap_bp(spans, starts, chrom, s, e):
    lst = spans.get(chrom)
    if not lst:
        return 0
    i = max(bisect.bisect_right(starts[chrom], s) - 1, 0)
    total = 0
    while i < len(lst) and lst[i][0] <= e:
        total += max(0, min(e, lst[i][1]) - max(s, lst[i][0]) + 1)
        i += 1
    return total


def te_fraction_by_gene(genes, mrna_gene, cds, te_path):
    """Fraction of the CDS of the main isoform (longest CDS) lying in a TE."""
    spans, starts = read_te_spans(te_path)
    main = {}
    for mid, gid in mrna_gene.items():
        length = sum(e - s + 1 for _, s, e, _ in cds.get(mid, ()))
        if gid not in main or length > main[gid][0]:
            main[gid] = (length, mid)
    frac = {}
    for gid, (length, mid) in main.items():
        if length > 0:
            ov = sum(overlap_bp(spans, starts, c, s, e) for c, s, e, _ in cds[mid])
            frac[gid] = ov / length
    return frac


def reference_overlap_genes(genes, mrna_gene, cds, ref_gff3, min_frac=0.5):
    """Genes with >= min_frac of their main-isoform CDS overlapping a CDS of ref_gff3 on the same strand."""
    with open(ref_gff3) as fh:
        raw = defaultdict(list)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not line.startswith("#") and len(f) > 8 and f[2] == "CDS":
                raw[(f[0], f[6])].append((int(f[3]), int(f[4])))
    spans, starts = {}, {}
    for k, lst in raw.items():
        lst.sort()
        merged = []
        for s_, e_ in lst:
            if merged and s_ <= merged[-1][1] + 1:
                merged[-1][1] = max(merged[-1][1], e_)
            else:
                merged.append([s_, e_])
        spans[k] = merged
        starts[k] = [m[0] for m in merged]
    main = {}
    for mid, gid in mrna_gene.items():
        length = sum(e - s + 1 for _, s, e, _ in cds.get(mid, ()))
        if gid not in main or length > main[gid][0]:
            main[gid] = (length, mid)
    out = set()
    for gid, (length, mid) in main.items():
        if length:
            ov = sum(overlap_bp(spans, starts, (c, st), s, e) for c, s, e, st in cds[mid])
            if ov / length >= min_frac:
                out.add(gid)
    return out


def liftoff_carried(path):
    with open(path, newline="") as fh:
        return {r["final_gene_id"] for r in csv.DictReader(fh, delimiter="\t")
                if r.get("decision") == "carried_over_from_liftoff"}


def busco_genes(path):
    out = set()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not line.startswith("#") and len(f) >= 3 and f[2].strip():
                out.add(gene_of(f[2].strip()))
    return out


def host_function_genes(d2go, eggnog, ipr):
    """Genes with a functional hit that is not only a TE-encoded protein."""
    strict, te = set(), set()
    with open(d2go) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            g = gene_of(f[0])
            if len(f) > 4 and f[4] not in ("", "-"):
                strict.add(g)
            if len(f) > 1 and TE_TEXT.search(f[1]):
                te.add(g)
    with open(eggnog) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            g = gene_of(f[0])
            desc = f[7] if len(f) > 7 else "-"
            if (desc not in ("", "-") or (len(f) > 9 and f[9] not in ("", "-"))
                    or (len(f) > 11 and f[11] not in ("", "-"))):
                strict.add(g)
            if TE_TEXT.search(desc):
                te.add(g)
    with open(ipr) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            g = gene_of(f[0])
            if f[3] not in NON_INFORMATIVE_IPR:
                strict.add(g)
            text = (f[5] if len(f) > 5 else "") + " | " + (f[12] if len(f) > 12 else "")
            if TE_TEXT.search(text):
                te.add(g)
    return strict - te


def filter_gff3(src, dst, drop):
    """Write every line except those hanging (at any depth) below a dropped gene."""
    parents = {}
    with open(src) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 9:
                a = parse_attrs(f[8])
                if "ID" in a and "Parent" in a:
                    parents[a["ID"]] = a["Parent"].split(",")
    memo = {}

    def below(fid):
        if fid in drop:
            return True
        if fid not in memo:
            memo[fid] = False
            memo[fid] = any(below(p) for p in parents.get(fid, ()))
        return memo[fid]

    kept = 0
    with open(src) as fh, open(dst, "w") as out:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(f) < 9:
                out.write(line)
                continue
            a = parse_attrs(f[8])
            ps = a.get("Parent", "").split(",") if "Parent" in a else []
            if below(a.get("ID", "")) or (ps and all(below(p) for p in ps)):
                continue
            kept += f[2] == "gene"
            out.write(line)
    return kept


def filter_fasta(src, dst, drop):
    n = 0
    keep = True
    with open(src) as fin, open(dst, "w") as out:
        for line in fin:
            if line.startswith(">"):
                keep = gene_of(line[1:].split()[0]) not in drop
                n += keep
            if keep:
                out.write(line)
    return n


def filter_table(src, dst, drop):
    """Drop rows whose first column (protein/CDS-record or gene ID) belongs to a dropped gene.
    '#' lines (headers/comments) are always kept."""
    with open(src) as fin, open(dst, "w") as out:
        for line in fin:
            if line.startswith("#") or gene_of(line.split("\t", 1)[0]) not in drop:
                out.write(line)


MIN_HOMOL_PIDENT = 35.0   # same thresholds as scripts/classify_unsupported_genes.py
MIN_HOMOL_QCOV = 0.5
DIAMOND_COLS = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore".split()


def select_candidates(a):
    genes, mrna_gene, cds = read_gff3(a.gff3)
    te_frac = te_fraction_by_gene(genes, mrna_gene, cds, a.te_gff3)
    protected = liftoff_carried(a.liftoff_correspondence_tsv) | busco_genes(a.busco_full_table)
    protected |= host_function_genes(a.diamond2go_main, a.eggnog_main, a.interproscan_main)
    protected |= reference_overlap_genes(genes, mrna_gene, cds, a.reference_gff3)
    cand = {g for g, fr in te_frac.items() if fr >= a.te_frac and g not in protected}
    with open(os.path.join(a.out_dir, "te_candidates.tsv"), "w") as fh:
        fh.write("gene_id\tte_frac_cds\n")
        for g in sorted(cand):
            fh.write("%s\t%.3f\n" % (g, te_frac[g]))
    # candidate main proteins, gene-ID headers, '*' / '-' stripped (DIAMOND / HMMER-safe)
    n = 0
    keep = False
    with open(a.proteins_main) as fin, open(os.path.join(a.out_dir, "te_candidates.fasta"), "w") as out:
        for line in fin:
            if line.startswith(">"):
                g = gene_of(line[1:].split()[0])
                keep = g in cand
                if keep:
                    out.write(">%s\n" % g)
                    n += 1
            elif keep:
                out.write(re.sub(r"[*-]", "", line))
    print("%d TE-derived candidate genes (%d proteins) of %d genes" % (len(cand), n, len(genes)))


TE_FUNCTION_SOURCES = ("InterPro", "Pfam", "PANTHER", "eggNOG", "Diamond2GO")


def read_te_function_terms(path):
    """source -> set of terms; TSV with columns source, term, name (term = accession for InterPro / Pfam /
    PANTHER, the exact description in 'name' for eggNOG / Diamond2GO)."""
    terms = {s: set() for s in TE_FUNCTION_SOURCES}
    with open(path, newline="") as fh:
        for r in csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"):
            src = r["source"].strip()
            if src not in terms:
                raise SystemExit("unknown source %r in %s (expected one of %s)" % (src, path, ", ".join(TE_FUNCTION_SOURCES)))
            terms[src].add(r["term"].strip() if src in ("InterPro", "Pfam", "PANTHER") else r["name"].strip())
    return terms


def te_function_genes(terms, d2go, eggnog, ipr):
    """Genes whose main protein carries at least one of the listed TE functional annotations."""
    out = set()
    with open(ipr) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            if (len(f) > 11 and f[11] in terms["InterPro"]) or (f[3] in ("Pfam", "PANTHER") and f[4] in terms[f[3]]):
                out.add(gene_of(f[0]))
    with open(eggnog) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not line.startswith("#") and len(f) > 7 and f[7] in terms["eggNOG"]:
                out.add(gene_of(f[0]))
    with open(d2go) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("#") or len(f) < 2:
                continue
            desc = f[1].strip()
            desc = desc.split(" ", 1)[1].strip() if " " in desc else desc   # drop the hit accession
            if desc in terms["Diamond2GO"]:
                out.add(gene_of(f[0]))
    return out


def nonvitis_homolog_genes(hits_path):
    best = {}
    with open(hits_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(DIAMOND_COLS):
                continue
            q, pident, qlen, qs, qe, bits = f[0], float(f[2]), int(f[4]), int(f[6]), int(f[7]), float(f[11])
            if q not in best or bits > best[q][0]:
                best[q] = (bits, pident, (qe - qs + 1) / qlen)
    return {q for q, (_, pid, qcov) in best.items() if pid >= MIN_HOMOL_PIDENT and qcov >= MIN_HOMOL_QCOV}


def finalize(a):
    out = a.out_dir
    fdir = os.path.join(out, "functional_annotation")
    for sub in ("diamond2go", "eggnog", "interproscan"):
        os.makedirs(os.path.join(fdir, sub), exist_ok=True)

    genes, mrna_gene, _cds = read_gff3(a.gff3)
    coding = set(mrna_gene.values())
    ncrna = set(genes) - coding

    cand = {}
    with open(a.candidates_tsv) as fh:
        next(fh)
        for line in fh:
            g, fr = line.rstrip("\n").split("\t")
            cand[g] = float(fr)
    rescued = nonvitis_homolog_genes(a.diamond_hits) & set(cand)
    te_drop = set(cand) - rescued
    fn_drop, n_terms = set(), 0
    if a.te_function_terms:
        if not a.liftoff_correspondence_tsv:
            raise SystemExit("--te-function-terms needs --liftoff-correspondence-tsv (Liftoff genes are kept)")
        terms = read_te_function_terms(a.te_function_terms)
        n_terms = sum(len(v) for v in terms.values())
        fn_drop = (te_function_genes(terms, a.diamond2go_main, a.eggnog_main, a.interproscan_main) & coding) \
            - liftoff_carried(a.liftoff_correspondence_tsv) - te_drop
    drop = te_drop | fn_drop | ncrna

    with open(os.path.join(out, "genes_removed.tsv"), "w") as fh:
        fh.write("gene_id\treason\tte_frac_cds\n")
        for g in sorted(drop):
            reason = "te_derived" if g in te_drop else "te_function" if g in fn_drop else "ncRNA"
            fh.write("%s\t%s\t%s\n" % (g, reason,
                                       "%.3f" % cand[g] if g in cand else ""))

    n_genes = filter_gff3(a.gff3, os.path.join(out, "final_annotation.noTE.gff3"), drop)
    n_main = filter_fasta(a.proteins_main, os.path.join(out, "final_annotation_proteins_main.noTE.fasta"), drop)
    n_all = filter_fasta(a.proteins_all, os.path.join(out, "final_annotation_proteins_all.noTE.fasta"), drop)
    for tag, d2, eg, ip in (("all", a.diamond2go_all, a.eggnog_all, a.interproscan_all),
                            ("main", a.diamond2go_main, a.eggnog_main, a.interproscan_main)):
        filter_table(d2, os.path.join(fdir, "diamond2go", "final_annotation_proteins_%s.noTE.diamond2go.tsv" % tag), drop)
        filter_table(eg, os.path.join(fdir, "eggnog", "final_annotation_proteins_%s.noTE.emapper.annotations" % tag), drop)
        filter_table(ip, os.path.join(fdir, "interproscan", "final_annotation_proteins_%s.noTE.tsv" % tag), drop)

    summary = {
        "schema_version": "titan.primary_noTE.v2",
        "te_frac_cds_threshold": a.te_frac,
        "genes_in": len(genes), "coding_genes_in": len(coding),
        "te_candidates": len(cand),
        "te_candidates_rescued_by_nonvitis_homolog": len(rescued),
        "te_genes_removed": len(te_drop),
        "te_function_terms": n_terms,
        "te_function_genes_removed": len(fn_drop),
        "ncrna_genes_removed": len(ncrna),
        "genes_out": n_genes,
        "proteins_main_out": n_main, "proteins_all_out": n_all,
    }
    with open(os.path.join(out, "primary_noTE_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(json.dumps(summary, sort_keys=True))


def candidate_structure(gff3):
    """Per gene: (type, number of exons of the main isoform, CDS length of the main isoform)."""
    genes, mrna_type, mrna_gene = set(), {}, {}
    exons, cds_len = defaultdict(int), defaultdict(int)
    with open(gff3) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = parse_attrs(f[8])
            if f[2] == "gene":
                genes.add(a["ID"])
            elif f[2] in ("mRNA", "ncRNA", "transcript") and "Parent" in a:
                mrna_type[a["ID"]] = f[2]
                mrna_gene[a["ID"]] = a["Parent"].split(",")[0]
            elif f[2] == "exon" and "Parent" in a:
                for q in a["Parent"].split(","):
                    exons[q] += 1
            elif f[2] == "CDS" and "Parent" in a:
                for q in a["Parent"].split(","):
                    cds_len[q] += int(f[4]) - int(f[3]) + 1
    main = {}
    for mid, g in mrna_gene.items():
        if mrna_type[mid] != "mRNA":
            continue
        if g not in main or cds_len[mid] > main[g][2]:
            main[g] = (mid, exons[mid], cds_len[mid])
    ncrna = {g for mid, g in mrna_gene.items() if mrna_type[mid] == "ncRNA"} - set(main)
    return genes, main, ncrna, sum(1 for t in mrna_type.values() if t == "mRNA")


def parse_busco_summary(path):
    """C/S/D/F/M percentages of a BUSCO short summary ('' when absent or empty)."""
    keys = ("C", "S", "D", "F", "M")
    if path in (None, "NA") or not os.path.exists(path):
        return dict.fromkeys(keys, "")
    for line in open(path):
        m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%", line)
        if m:
            return dict(zip(keys, m.groups()))
    return dict.fromkeys(keys, "")


def parse_omark_sum(path):
    """Completeness (single/duplicated/missing HOGs) and consistency (A/I/C/U) percentages of an OMArk .sum."""
    keys = ("omark_S", "omark_D", "omark_M", "omark_consistent", "omark_inconsistent",
            "omark_contaminant", "omark_unknown")
    out = dict.fromkeys(keys, "")
    if path in (None, "NA") or not os.path.exists(path):
        return out
    pct = [l for l in open(path) if not l.startswith("#") and "%" in l]
    if len(pct) >= 2:
        n = lambda tag, line: (re.search(tag + r":([\d.]+)%", line) or [None, ""])[1]
        out.update(omark_S=n("S", pct[0]), omark_D=n("D", pct[0]), omark_M=n("M", pct[0]),
                   omark_consistent=n("A", pct[1]), omark_inconsistent=n("I", pct[1]),
                   omark_contaminant=n("C", pct[1]), omark_unknown=n("U", pct[1]))
    return out


def compare(a):
    unsupported = set(json.load(open(a.expression_json)).get("unsupported_gene_ids", []))
    with_function = set()
    for path in (a.diamond2go_main, a.eggnog_main, a.interproscan_main):
        with open(path) as fh:
            for line in fh:
                if line.strip() and not line.startswith("#"):
                    with_function.add(gene_of(line.split("\t", 1)[0]))
    spans, starts = read_te_spans(a.te_gff3)
    rows = []
    for label, gff3, busco, omark in a.candidate:
        genes, main, ncrna, n_mrna = candidate_structure(gff3)
        _g, mrna_gene, cds = read_gff3(gff3)
        te = te_fraction_by_gene(_g, mrna_gene, cds, a.te_gff3)
        n = len(main) or 1
        lens = sorted(v[2] // 3 for v in main.values())
        pct = lambda k: "%.1f" % (100.0 * k / n)
        row = {
            "candidate": label, "genes": len(genes), "protein_coding_genes": len(main),
            "ncrna_genes": len(ncrna), "mrna": n_mrna,
            "mono_exonic_pct": pct(sum(1 for v in main.values() if v[1] <= 1)),
            "median_protein_aa": lens[len(lens) // 2] if lens else "",
            "te_cds_ge50_pct": pct(sum(1 for g in main if te.get(g, 0) >= 0.5)),
            "te_cds_ge80_pct": pct(sum(1 for g in main if te.get(g, 0) >= 0.8)),
            "with_function_pct": pct(sum(1 for g in main if g in with_function)),
            "expressed_pct": pct(sum(1 for g in main if g not in unsupported)),
        }
        row.update({"busco_" + k: v for k, v in parse_busco_summary(busco).items()})
        row.update(parse_omark_sum(omark))
        rows.append(row)
    cols = list(rows[0])
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(open(a.out).read())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    sub.required = True
    s = sub.add_parser("select", help="list TE-derived candidate genes and write their proteins")
    f = sub.add_parser("finalize", help="apply the filters and write the primary_noTE files")
    for p in (s, f):
        p.add_argument("--gff3", required=True)
        p.add_argument("--proteins-main", required=True)
        p.add_argument("--te-frac", type=float, default=0.80)
        p.add_argument("--out-dir", default=".")
    s.add_argument("--te-gff3", required=True)
    s.add_argument("--reference-gff3", required=True, help="Liftoff-transferred previous annotation")
    s.add_argument("--liftoff-correspondence-tsv", required=True)
    s.add_argument("--busco-full-table", required=True)
    f.add_argument("--candidates-tsv", required=True)
    f.add_argument("--diamond-hits", required=True, help="DIAMOND outfmt 6 (%s)" % " ".join(DIAMOND_COLS))
    f.add_argument("--proteins-all", required=True)
    f.add_argument("--te-function-terms", help="TSV (source, term, name) of TE functional annotations: genes "
                   "carrying one are removed unless carried over from Liftoff")
    f.add_argument("--liftoff-correspondence-tsv", dest="liftoff_correspondence_tsv",
                   help="AEGIS liftoff_gene_id_correspondence.tsv (required with --te-function-terms)")
    for p in (s, f):
        p.add_argument("--diamond2go-main", required=True)
        p.add_argument("--eggnog-main", required=True, help="*.emapper.annotations")
        p.add_argument("--interproscan-main", required=True, help="InterProScan TSV")
    f.add_argument("--diamond2go-all", required=True)
    f.add_argument("--eggnog-all", required=True)
    f.add_argument("--interproscan-all", required=True)
    c = sub.add_parser("compare", help="side-by-side quality table of candidate annotations")
    c.add_argument("--candidate", nargs=4, action="append", required=True,
                   metavar=("LABEL", "GFF3", "BUSCO_SHORT_SUMMARY", "OMARK_SUM"),
                   help="repeatable; give NA for a missing BUSCO/OMArk file")
    c.add_argument("--te-gff3", required=True)
    c.add_argument("--expression-json", required=True, help="expression_support_summary.json")
    c.add_argument("--diamond2go-main", required=True)
    c.add_argument("--eggnog-main", required=True)
    c.add_argument("--interproscan-main", required=True)
    c.add_argument("--out", default="candidates_comparison.tsv")
    a = ap.parse_args()
    if a.cmd == "compare":
        compare(a)
        return
    os.makedirs(a.out_dir, exist_ok=True)
    select_candidates(a) if a.cmd == "select" else finalize(a)


if __name__ == "__main__":
    main()
