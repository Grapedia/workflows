#!/usr/bin/env python3
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "soft_mask_genome.py"

ORIGINAL = ">chr1 desc\nACGTACGTAC\nGTACGT\n>chr2\nNNACGTAC\n"
# EDTA order differs (chr2 first), repeats at chr1[2:6] and chr2[4:6]; chr2 starts with a native gap
MASKED = ">chr2\nNNACNNAC\n>chr1\nACNNNNGTAC\nGTACGT\n"


def run(original, masked):
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "o.fa").write_text(original)
        (tmp / "m.fa").write_text(masked)
        proc = subprocess.run(
            [sys.executable, str(SCRIPT), "--genome", str(tmp / "o.fa"),
             "--masked", str(tmp / "m.fa"), "--out", str(tmp / "s.fa")],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        out = (tmp / "s.fa").read_text() if (tmp / "s.fa").exists() else None
        return proc, out


def main():
    proc, out = run(ORIGINAL, MASKED)
    assert proc.returncode == 0, proc.stderr
    assert out == ">chr1 desc\nACgtacGTACGTACGT\n>chr2\nNNACgtAC\n", repr(out)

    proc, _ = run(ORIGINAL, MASKED.replace("GTACGT\n", "GTAC\n"))
    assert proc.returncode != 0 and "length differs" in proc.stderr, proc.stderr

    proc, _ = run(ORIGINAL, ">chr2\nNNACNNAC\n")
    assert proc.returncode != 0 and "not in" in proc.stderr, proc.stderr

    proc, _ = run(ORIGINAL, MASKED + ">chr3\nAC\n")
    assert proc.returncode != 0 and "only in" in proc.stderr, proc.stderr
    print("ok")


if __name__ == "__main__":
    main()
