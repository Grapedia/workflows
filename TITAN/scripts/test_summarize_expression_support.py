#!/usr/bin/env python3
import json
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "summarize_expression_support.py"


def test_expression_support_summary():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        gff = tmp / "annotation.gff3"
        gff.write_text(
            "\n".join(
                [
                    "##gff-version 3",
                    "chr1\tTITAN\tgene\t1\t100\t.\t+\t.\tID=gene1",
                    "chr1\tTITAN\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1",
                    "chr1\tTITAN\tgene\t200\t300\t.\t+\t.\tID=gene2",
                    "chr1\tTITAN\tmRNA\t200\t300\t.\t+\t.\tID=tx2;Parent=gene2",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        quant = tmp / "sampleA_quant"
        quant.mkdir()
        (quant / "quant.sf").write_text(
            "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
            "tx1\t100\t90\t2.5\t25\n"
            "tx2\t100\t90\t0.1\t1\n",
            encoding="utf-8",
        )

        output_json = tmp / "summary.json"
        output_mqc = tmp / "summary_mqc.tsv"
        output_matrix = tmp / "matrix.tsv"
        subprocess.run(
            [
                sys.executable,
                str(SCRIPT),
                "--gff",
                str(gff),
                "--quant-dirs",
                str(quant),
                "--min-tpm",
                "0.5",
                "-o",
                str(output_json),
                "--multiqc-tsv",
                str(output_mqc),
                "--tpm-matrix",
                str(output_matrix),
            ],
            check=True,
        )

        summary = json.loads(output_json.read_text(encoding="utf-8"))
        assert summary["total_genes"] == 2
        assert summary["supported_genes"] == 1
        assert summary["unsupported_gene_ids"] == ["gene2"]
        assert "Supported genes (%)\t50.0" in output_mqc.read_text(encoding="utf-8")


if __name__ == "__main__":
    test_expression_support_summary()
