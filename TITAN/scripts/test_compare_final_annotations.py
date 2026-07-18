#!/usr/bin/env python3
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "compare_final_annotations.py"


def test_compare_final_annotations_counts_overlaps():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        aegis = tmp / "aegis.gff3"
        mikado = tmp / "mikado.gff3"
        aegis.write_text("##gff-version 3\nchr1\tA\tgene\t10\t50\t.\t+\t.\tID=a1\n", encoding="utf-8")
        mikado.write_text(
            "##gff-version 3\n"
            "chr1\tM\tgene\t30\t70\t.\t+\t.\tID=m1\n"
            "chr1\tM\tgene\t100\t130\t.\t+\t.\tID=m2\n",
            encoding="utf-8",
        )
        subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "--aegis-gff3",
                str(aegis),
                "--mikado-gff3",
                str(mikado),
            ],
            cwd=tmp,
            check=True,
        )
        report = (tmp / "final_annotation_sources.json").read_text(encoding="utf-8")
        assert '"mikado_gene_count": 2' in report
        assert '"mikado_genes_overlapping_aegis": 1' in report


if __name__ == "__main__":
    test_compare_final_annotations_counts_overlaps()
    print("compare_final_annotations tests OK")
