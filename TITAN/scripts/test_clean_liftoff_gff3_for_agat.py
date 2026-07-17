#!/usr/bin/env python3
"""Tests for clean_liftoff_gff3_for_agat.py."""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "clean_liftoff_gff3_for_agat.py"


def run_cleaner(input_text: str) -> tuple[str, str]:
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        input_gff3 = tmpdir_path / "input.gff3"
        output_gff3 = tmpdir_path / "cleaned.gff3"
        removed_ids = tmpdir_path / "removed_ids.txt"
        input_gff3.write_text(input_text, encoding="utf-8")

        subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "--input",
                str(input_gff3),
                "--output",
                str(output_gff3),
                "--removed-ids",
                str(removed_ids),
            ],
            check=True,
        )

        return output_gff3.read_text(encoding="utf-8"), removed_ids.read_text(encoding="utf-8")


def test_deleted_obsolete_feature_and_children_are_removed() -> None:
    cleaned, removed = run_cleaner(
        "\n".join(
            [
                "##gff-version 3",
                "chr1\tLiftoff\tgene\t1\t100\t.\t+\t.\tID=gene_deleted;obsolete=true;Note=Deleted in source",
                "chr1\tLiftoff\tmRNA\t1\t100\t.\t+\t.\tID=mrna_deleted;Parent=gene_deleted",
                "chr1\tLiftoff\tCDS\t1\t50\t.\t+\t0\tID=cds_deleted;Parent=mrna_deleted",
                "chr1\tLiftoff\tgene\t200\t300\t.\t+\t.\tID=gene_kept",
                "chr1\tLiftoff\tmRNA\t200\t300\t.\t+\t.\tID=mrna_kept;Parent=gene_kept",
                "malformed line that AGAT should not receive",
                "",
            ]
        )
    )

    assert "gene_deleted" not in cleaned
    assert "mrna_deleted" not in cleaned
    assert "cds_deleted" not in cleaned
    assert "gene_kept" in cleaned
    assert "mrna_kept" in cleaned
    assert "malformed line" not in cleaned
    assert removed.splitlines() == ["cds_deleted", "gene_deleted", "mrna_deleted"]


def test_empty_removal_report_is_created_when_no_feature_matches() -> None:
    cleaned, removed = run_cleaner(
        "\n".join(
            [
                "##gff-version 3",
                "chr1\tLiftoff\tgene\t1\t100\t.\t+\t.\tID=gene_kept",
                "chr1\tLiftoff\tmRNA\t1\t100\t.\t+\t.\tID=mrna_kept;Parent=gene_kept",
                "",
            ]
        )
    )

    assert "gene_kept" in cleaned
    assert "mrna_kept" in cleaned
    assert removed == ""


def main() -> None:
    test_deleted_obsolete_feature_and_children_are_removed()
    test_empty_removal_report_is_created_when_no_feature_matches()
    print("Clean Liftoff GFF3 tests OK")


if __name__ == "__main__":
    main()
