#!/usr/bin/env python3

import tempfile
import textwrap
import unittest
from pathlib import Path

from build_egapx_paramfile import build_paramfile


class BuildEgapxParamfileTests(unittest.TestCase):
    def test_replaces_only_short_reads_with_trimmed_paths(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            source = tmp / "input.yaml"
            manifest = tmp / "trimmed.tsv"
            output = tmp / "derived.yaml"
            read_1 = tmp / "sample_1.trimmed.fastq.gz"
            read_2 = tmp / "sample_2.trimmed.fastq.gz"

            source.write_text(
                textwrap.dedent(
                    """\
                    genome: target.fa
                    taxid: 0000
                    organism: Synthetic test organism
                    short_reads:
                      - - stale_sample
                        - - stale_1.fastq.gz
                          - stale_2.fastq.gz
                    long_reads:
                      - - isoforms
                        - - hq_transcripts.fasta
                    """
                ),
                encoding="utf-8",
            )
            manifest.write_text(f"sample\tpaired\t{read_1}\t{read_2}\n", encoding="utf-8")

            build_paramfile(source, manifest, output)

            derived = output.read_text(encoding="utf-8")
            self.assertIn("genome: target.fa", derived)
            self.assertIn("long_reads:", derived)
            self.assertNotIn("stale_sample", derived)
            self.assertIn("sample_1.trimmed.fastq.gz", derived)
            self.assertIn("sample_2.trimmed.fastq.gz", derived)

    def test_single_end_rows_emit_one_read(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            source = tmp / "input.yaml"
            manifest = tmp / "trimmed.tsv"
            output = tmp / "derived.yaml"
            read_1 = tmp / "sample_1.trimmed.fastq.gz"

            source.write_text("genome: target.fa\ntaxid: 0000\norganism: Synthetic\n", encoding="utf-8")
            manifest.write_text(f"sample\tsingle\t{read_1}\t\n", encoding="utf-8")

            build_paramfile(source, manifest, output)

            derived = output.read_text(encoding="utf-8")
            self.assertIn("    - - " + str(read_1.resolve()), derived)
            self.assertNotIn("sample_2.trimmed.fastq.gz", derived)


if __name__ == "__main__":
    unittest.main()
