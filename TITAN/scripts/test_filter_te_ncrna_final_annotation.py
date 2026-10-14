#!/usr/bin/env python3
"""End-to-end test of scripts/filter_te_ncrna_final_annotation.py (select + finalize) on a toy annotation."""
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "filter_te_ncrna_final_annotation.py"


def gene(gid, start, end, strand="+", kind="mRNA"):
    return [
        "chr1\tT\tgene\t%d\t%d\t.\t%s\t.\tID=%s" % (start, end, strand, gid),
        "chr1\tT\t%s\t%d\t%d\t.\t%s\t.\tID=%s_t001;Parent=%s" % (kind, start, end, strand, gid, gid),
        "chr1\tT\texon\t%d\t%d\t.\t%s\t.\tID=%s_e1;Parent=%s_t001" % (start, end, strand, gid, gid),
    ] + (["chr1\tT\tCDS\t%d\t%d\t.\t%s\t0\tID=%s_t001_CDS1;Parent=%s_t001" % (start, end, strand, gid, gid)]
         if kind == "mRNA" else [])


# g_te: in a TE, unprotected -> removed; g_te_homol: in a TE, non-Vitis homolog -> kept;
# g_te_lift: in a TE, carried from Liftoff -> kept; g_te_fn: in a TE, host function -> kept;
# g_ok: not in a TE -> kept; g_nc: ncRNA -> removed
GFF = "\n".join(["##gff-version 3"]
                + gene("g_te", 1000, 1099) + gene("g_te_homol", 2000, 2099) + gene("g_te_lift", 3000, 3099)
                + gene("g_te_fn", 4000, 4099) + gene("g_ok", 5000, 5099) + gene("g_nc", 6000, 6099, kind="ncRNA")) + "\n"
TE = "\n".join("chr1\tEDTA\tLTR_retrotransposon\t%d\t%d\t.\t+\t.\tID=te%d" % (s, s + 99, s) for s in (1000, 2000, 3000, 4000)) + "\n"
REF = "chr1\tL\tCDS\t9000\t9100\t.\t+\t0\tID=x;Parent=y\n"
LIFT = "original_gene_id\tfinal_gene_id\tdecision\ng_te_lift\tg_te_lift\tcarried_over_from_liftoff\ng_te\tg_te\tnew_vitvi_id\n"
BUSCO = "# Busco id\tStatus\tSequence\n1at1\tMissing\t\n"
D2GO = "#gene_id\tgene_name\tspecies\te-value\tGO-term\n"
EGG = "#query\tseed\tevalue\tscore\tOGs\tmax_lvl\tCOG\tDescription\tPreferred_name\tGOs\n" \
      "g_te_fn_t001_CDS1.prot\ts\t0\t1\tog\tl\tK\tkinase activity\t-\t-\n"
IPR = "g_te_fn_t001_CDS1.prot\tmd5\t33\tMobiDBLite\tx\tdisorder\t1\t9\t-\tT\td\t-\t-\n"
HITS = "g_te_homol\tsp|A|B\t80.0\t30\t33\t300\t1\t33\t1\t33\t1e-9\t60\n"


def run(*args):
    subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], check=True, stdout=subprocess.DEVNULL)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        t = Path(tmp)
        files = {"a.gff3": GFF, "te.gff3": TE, "ref.gff3": REF, "lift.tsv": LIFT, "busco.tsv": BUSCO,
                 "d2go.tsv": D2GO, "egg.tsv": EGG, "ipr.tsv": IPR, "hits.tsv": HITS}
        for name, text in files.items():
            (t / name).write_text(text)
        prot = "".join(">%s_t001_CDS1.prot\nMAAA*\n" % g for g in ("g_te", "g_te_homol", "g_te_lift", "g_te_fn", "g_ok"))
        (t / "prot.fa").write_text(prot)
        common = ["--gff3", t / "a.gff3", "--proteins-main", t / "prot.fa", "--out-dir", t / "out",
                  "--diamond2go-main", t / "d2go.tsv", "--eggnog-main", t / "egg.tsv",
                  "--interproscan-main", t / "ipr.tsv"]
        run("select", *common, "--te-gff3", t / "te.gff3", "--reference-gff3", t / "ref.gff3",
            "--liftoff-correspondence-tsv", t / "lift.tsv", "--busco-full-table", t / "busco.tsv")
        cands = [l.split("\t")[0] for l in (t / "out" / "te_candidates.tsv").read_text().splitlines()[1:]]
        assert sorted(cands) == ["g_te", "g_te_homol"], cands
        assert (t / "out" / "te_candidates.fasta").read_text() == ">g_te\nMAAA\n>g_te_homol\nMAAA\n"
        run("finalize", *common, "--candidates-tsv", t / "out" / "te_candidates.tsv", "--diamond-hits", t / "hits.tsv",
            "--proteins-all", t / "prot.fa", "--diamond2go-all", t / "d2go.tsv", "--eggnog-all", t / "egg.tsv",
            "--interproscan-all", t / "ipr.tsv")
        out = (t / "out" / "final_annotation.noTE.gff3").read_text()
        genes = sorted(l.split("ID=")[1] for l in out.splitlines() if "\tgene\t" in l)
        assert genes == ["g_ok", "g_te_fn", "g_te_homol", "g_te_lift"], genes
        assert "g_nc" not in out and "g_te_t001" not in out
        removed = dict(l.split("\t")[:2] for l in (t / "out" / "genes_removed.tsv").read_text().splitlines()[1:])
        assert removed == {"g_te": "te_derived", "g_nc": "ncRNA"}, removed
        fa = (t / "out" / "final_annotation_proteins_main.noTE.fasta").read_text()
        assert ">g_te_t001_CDS1.prot" not in fa and ">g_ok_t001_CDS1.prot" in fa
        assert (t / "out" / "functional_annotation" / "eggnog" / "final_annotation_proteins_main.noTE.emapper.annotations").exists()
        # compare: primary vs primary_noTE
        (t / "expr.json").write_text('{"unsupported_gene_ids": ["g_te", "g_ok"]}')
        (t / "busco.txt").write_text("\tC:99.0%[S:98.0%,D:1.0%],F:0.5%,M:0.5%,n:10\n")
        run("compare", "--candidate", "primary", t / "a.gff3", t / "busco.txt", "NA",
            "--candidate", "primary_noTE", t / "out" / "final_annotation.noTE.gff3", "NA", "NA",
            "--te-gff3", t / "te.gff3", "--expression-json", t / "expr.json", "--diamond2go-main", t / "d2go.tsv",
            "--eggnog-main", t / "egg.tsv", "--interproscan-main", t / "ipr.tsv", "--out", t / "cmp.tsv")
        rows = [l.split("\t") for l in (t / "cmp.tsv").read_text().splitlines()]
        table = {r[0]: dict(zip(rows[0], r)) for r in rows[1:]}
        assert table["primary"]["genes"] == "6" and table["primary"]["ncrna_genes"] == "1", table["primary"]
        assert table["primary_noTE"]["genes"] == "4" and table["primary_noTE"]["ncrna_genes"] == "0"
        assert table["primary"]["te_cds_ge80_pct"] == "80.0" and table["primary_noTE"]["te_cds_ge80_pct"] == "75.0"
        assert table["primary"]["busco_C"] == "99.0" and table["primary_noTE"]["busco_C"] == ""
    print("OK")


if __name__ == "__main__":
    main()
