#!/usr/bin/env python3
import io
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "trnascan_to_gff3.py"
FIXTURE = ROOT / "test-data" / "minimal" / "valid" / "trnascan.out"


spec = importlib.util.spec_from_file_location("trnascan_to_gff3", SCRIPT)
trnascan_to_gff3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(trnascan_to_gff3)


def test_fixture_conversion():
    output = io.StringIO()
    count = trnascan_to_gff3.convert(FIXTURE, output)
    lines = output.getvalue().strip().splitlines()
    features = [line for line in lines if not line.startswith("#")]

    assert count == 2
    assert lines[0] == "##gff-version 3"
    assert len(features) == 2
    assert all(len(feature.split("\t")) == 9 for feature in features)
    assert features[0].split("\t")[2] == "tRNA"
    assert features[0].split("\t")[6] == "+"
    assert features[1].split("\t")[3:7] == ["40", "110", "88.0", "-"]
    assert "product=tRNA-Gly" in features[1]


if __name__ == "__main__":
    test_fixture_conversion()
    print("trnascan_to_gff3 tests OK")
