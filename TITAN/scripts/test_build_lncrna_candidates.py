#!/usr/bin/env python3
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_lncrna_candidates.py"


def write(path, text):
    path.write_text(text, encoding="utf-8")


def test_candidate_builder_filters_overlaps_and_short_transcripts():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        genome = tmp / "genome.fa"
        final_gff3 = tmp / "final.gff3"
        trna = tmp / "trna.gff3"
        rfam = tmp / "rfam.gff3"
        gtf = tmp / "transcripts.gtf"

        write(genome, ">chr1\n" + "A" * 500 + "\n")
        write(final_gff3, "##gff-version 3\nchr1\tT\tCDS\t10\t50\t.\t+\t0\tID=cds1\n")
        write(trna, "##gff-version 3\nchr1\tT\ttRNA\t300\t330\t.\t+\t.\tID=t1\n")
        write(rfam, "##gff-version 3\n")
        write(
            gtf,
            "\n".join(
                [
                    'chr1\tT\ttranscript\t100\t220\t.\t+\t.\tgene_id "g1"; transcript_id "keep1";',
                    'chr1\tT\texon\t100\t220\t.\t+\t.\tgene_id "g1"; transcript_id "keep1";',
                    'chr1\tT\ttranscript\t1\t120\t.\t+\t.\tgene_id "g2"; transcript_id "coding_overlap";',
                    'chr1\tT\texon\t1\t120\t.\t+\t.\tgene_id "g2"; transcript_id "coding_overlap";',
                    'chr1\tT\ttranscript\t250\t360\t.\t+\t.\tgene_id "g3"; transcript_id "trna_overlap";',
                    'chr1\tT\texon\t250\t360\t.\t+\t.\tgene_id "g3"; transcript_id "trna_overlap";',
                    'chr1\tT\ttranscript\t400\t430\t.\t+\t.\tgene_id "g4"; transcript_id "short";',
                    'chr1\tT\texon\t400\t430\t.\t+\t.\tgene_id "g4"; transcript_id "short";',
                ]
            )
            + "\n",
        )

        subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "--genome",
                str(genome),
                "--final-annotation",
                str(final_gff3),
                "--trna-gff3",
                str(trna),
                "--rfam-gff3",
                str(rfam),
                "--min-length",
                "100",
                str(gtf),
            ],
            cwd=tmp,
            check=True,
        )
        gff3 = (tmp / "lncrna_candidates.gff3").read_text(encoding="utf-8")
        summary = (tmp / "lncrna_classification_summary.tsv").read_text(encoding="utf-8")
        assert "keep1" in gff3
        assert "coding_overlap" not in gff3
        assert "trna_overlap" not in gff3
        assert "lncRNA_candidate\t1" in summary
        assert "excluded_short\t1" in summary
        assert "excluded_overlap_coding_or_ncrna\t2" in summary
        assert "excluded_cpat_coding\t0" in summary

        cpat = tmp / "cpat_plant.output.ORF_prob.best.tsv"
        write(cpat, "seq_ID\tcoding_prob\nlncrna_candidate_1.t1\t0.73\n")
        subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "--genome",
                str(genome),
                "--final-annotation",
                str(final_gff3),
                "--trna-gff3",
                str(trna),
                "--rfam-gff3",
                str(rfam),
                "--min-length",
                "100",
                "--cpat-best-tsv",
                str(cpat),
                "--cpat-cutoff",
                "0.46",
                str(gtf),
            ],
            cwd=tmp,
            check=True,
        )
        gff3 = (tmp / "lncrna_candidates.gff3").read_text(encoding="utf-8")
        summary = (tmp / "lncrna_classification_summary.tsv").read_text(encoding="utf-8")
        assert "keep1" not in gff3
        assert "lncRNA_candidate\t0" in summary
        assert "excluded_cpat_coding\t1" in summary


if __name__ == "__main__":
    test_candidate_builder_filters_overlaps_and_short_transcripts()
    print("build_lncrna_candidates tests OK")
