#!/usr/bin/env python3
"""Smoke tests for shared shell helpers used by Nextflow modules."""

import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def assert_usage(script: str) -> None:
    path = ROOT / "scripts" / script
    result = subprocess.run([str(path)], cwd=ROOT, text=True, capture_output=True)
    assert result.returncode == 2, (script, result.returncode, result.stdout, result.stderr)
    assert "Usage:" in result.stderr, (script, result.stderr)


def main() -> None:
    assert_usage("run_aegis_merge.sh")
    assert_usage("run_stringtie_transcriptome.sh")
    print("Shared shell script tests OK")


if __name__ == "__main__":
    main()
