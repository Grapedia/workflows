#!/usr/bin/env python3
import io
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "rfam_tblout_to_gff3.py"
FIXTURE = ROOT / "test-data" / "minimal" / "valid" / "rfam_hits.tbl"


spec = importlib.util.spec_from_file_location("rfam_tblout_to_gff3", SCRIPT)
rfam_tblout_to_gff3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rfam_tblout_to_gff3)


def test_fixture_conversion():
    output = io.StringIO()
    count = rfam_tblout_to_gff3.convert(FIXTURE, output)
    lines = output.getvalue().strip().splitlines()
    features = [line for line in lines if not line.startswith("#")]

    assert count == 3
    assert lines[0] == "##gff-version 3"
    assert len(features) == 3
    assert all(len(feature.split("\t")) == 9 for feature in features)
    assert features[0].split("\t")[2] == "rRNA"
    assert features[1].split("\t")[2] == "snRNA"
    assert features[2].split("\t")[2] == "snoRNA"
    assert features[2].split("\t")[3:7] == ["40", "110", "55.5", "-"]
    assert "Rfam_ID=RF00012" in features[1]


if __name__ == "__main__":
    test_fixture_conversion()
    print("rfam_tblout_to_gff3 tests OK")
