#!/usr/bin/env python3
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "make_mikado_list.py"


def test_make_mikado_list_skips_empty_inputs_and_keeps_labels():
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        gtf = tmp / "input.gtf"
        empty = tmp / "empty.gtf"
        output = tmp / "list.tsv"
        gtf.write_text('chr1\tT\ttranscript\t1\t10\t.\t+\t.\tgene_id "g"; transcript_id "t";\n', encoding="utf-8")
        empty.write_text("", encoding="utf-8")

        subprocess.run(
            [
                "python3",
                str(SCRIPT),
                "--source",
                f"{gtf}:star_stringtie:True:5:False",
                "--source",
                f"{empty}:empty:False:1:False",
                "-o",
                str(output),
            ],
            check=True,
        )

        assert output.read_text(encoding="utf-8") == f"{gtf}\tstar_stringtie\tTrue\t5\tFalse\n"


if __name__ == "__main__":
    test_make_mikado_list_skips_empty_inputs_and_keeps_labels()
    print("make_mikado_list tests OK")
