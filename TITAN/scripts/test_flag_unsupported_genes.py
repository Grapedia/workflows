#!/usr/bin/env python3
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "flag_unsupported_genes.py"

GFF = "\n".join([
    "##gff-version 3",
    "chr1\tT\tgene\t1\t100\t.\t+\t.\tID=g_expr",
    "chr1\tT\tmRNA\t1\t100\t.\t+\t.\tID=g_expr_t001;Parent=g_expr",
    "chr1\tT\texon\t1\t100\t.\t+\t.\tID=g_expr_e001;Parent=g_expr_t001",
    "chr1\tT\tCDS\t1\t100\t.\t+\t0\tID=g_expr_t001_CDS1;Parent=g_expr_t001",
    "chr1\tT\tgene\t200\t300\t.\t+\t.\tID=g_nosig",
    "chr1\tT\tmRNA\t200\t300\t.\t+\t.\tID=g_nosig_t001;Parent=g_nosig",
    "chr1\tT\texon\t200\t300\t.\t+\t.\tID=g_nosig_e001;Parent=g_nosig_t001",
    "chr1\tT\tCDS\t200\t300\t.\t+\t0\tID=g_nosig_t001_CDS1;Parent=g_nosig_t001",
    "chr1\tT\tgene\t400\t500\t.\t-\t.\tID=g_lift",
    "chr1\tT\tmRNA\t400\t500\t.\t-\t.\tID=g_lift_t001;Parent=g_lift",
    "chr1\tT\texon\t400\t500\t.\t-\t.\tID=g_lift_e001;Parent=g_lift_t001",
    "chr1\tT\tgene\t600\t700\t.\t-\t.\tID=g_te",
    "chr1\tT\tmRNA\t600\t700\t.\t-\t.\tID=g_te_t001;Parent=g_te",
    "chr1\tT\texon\t600\t700\t.\t-\t.\tID=g_te_e001;Parent=g_te_t001",
    "chr1\tT\tgene\t800\t900\t.\t+\t.\tID=g_ncrna",
    "chr1\tT\tncRNA\t800\t900\t.\t+\t.\tID=g_ncrna_t001;Parent=g_ncrna",
    "chr1\tT\texon\t800\t900\t.\t+\t.\tID=g_ncrna_e001;Parent=g_ncrna_t001",
]) + "\n"
CLASSES = ("gene_id\tevidence_class\n"
           "g_expr\texpressed\ng_nosig\tno_signal\ng_lift\tliftoff_only\ng_te\tte_protein\n")
FASTA = ">g_expr_t001_CDS1.prot\nMKK\n>g_nosig_t001_CDS1.prot\nMAA\n>g_te_t001_CDS1.prot\nMTT\n"


def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT)] + [str(a) for a in args],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, check=True)


def genes(path):
    return sorted(l.split("\t")[8].split("ID=")[1].split(";")[0]
                  for l in Path(path).read_text().splitlines() if "\tgene\t" in l)


def test_mark_filter_policies():
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "in.gff3").write_text(GFF)
        (d / "classes.tsv").write_text(CLASSES)
        (d / "in.fa").write_text(FASTA)
        base = ["--gff3", d / "in.gff3", "--classes", d / "classes.tsv"]

        run("mark", *base, "--out", d / "marked.gff3")
        marked = (d / "marked.gff3").read_text()
        assert "ID=g_nosig;titan_evidence_class=no_signal;titan_filter_level=strict" in marked
        assert "ID=g_ncrna;titan_evidence_class=not_evaluated;titan_filter_level=keep" in marked
        assert marked.count("\n") == GFF.count("\n"), "mark must not add or drop lines"

        expected = {
            "strict": ["g_expr", "g_lift", "g_ncrna", "g_te"],
            "broad": ["g_expr", "g_ncrna", "g_te"],
            "broad+te": ["g_expr", "g_ncrna"],
        }
        for policy, kept in expected.items():
            out = d / ("f_%s.gff3" % policy.replace("+", "_"))
            fa = d / ("f_%s.fa" % policy.replace("+", "_"))
            run("filter", *base, "--policy", policy, "--out", out, "--proteins", d / "in.fa",
                "--proteins-out", fa)
            assert genes(out) == kept, (policy, genes(out))
            text = out.read_text()
            # no orphan child of a removed gene, header preserved
            assert text.startswith("##gff-version 3")
            for gid in ("g_expr", "g_nosig", "g_lift", "g_te", "g_ncrna"):
                has_gene = gid in kept
                assert (("ID=%s_t001" % gid) in text) == has_gene, (policy, gid)
                assert (("ID=%s_e001" % gid) in text) == has_gene, (policy, gid)
            fasta_ids = [l[1:].split("_t001")[0] for l in fa.read_text().splitlines() if l.startswith(">")]
            assert fasta_ids == [g for g in ("g_expr", "g_nosig", "g_te") if g in kept], (policy, fasta_ids)

        # the input is never modified
        assert (d / "in.gff3").read_text() == GFF

        out = run("summary", *base).stdout
        assert "strict" in out and "broad+te" in out


def test_unknown_class_is_rejected():
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "in.gff3").write_text(GFF)
        (d / "classes.tsv").write_text("gene_id\tevidence_class\ng_expr\tmystery\n")
        try:
            run("summary", "--gff3", d / "in.gff3", "--classes", d / "classes.tsv")
        except subprocess.CalledProcessError as exc:
            assert "unknown evidence class" in exc.stderr
        else:
            raise AssertionError("an unknown class must fail loudly")


if __name__ == "__main__":
    test_mark_filter_policies()
    test_unknown_class_is_rejected()
    print("ok")
