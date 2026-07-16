#!/usr/bin/env python3
import gzip
import hashlib
import io
import json
import sys
import tempfile
from contextlib import redirect_stderr
from pathlib import Path

import download_sra_fastq


def write_fastq_gz(path, name, sequence="ACGT"):
    with gzip.open(path, "wt") as handle:
        handle.write(f"@{name}\n{sequence}\n+\n{'!' * len(sequence)}\n")


def md5sum(path):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_ena_json(mock_dir, accession, layout, urls, md5s=None, sizes=None):
    record = {
        "run_accession": accession,
        "library_layout": layout,
        "fastq_ftp": ";".join(urls),
        "fastq_md5": ";".join(md5s or []),
        "fastq_bytes": ";".join(str(size) for size in (sizes or [])),
        "instrument_platform": "ILLUMINA",
        "library_strategy": "RNA-Seq",
    }
    with (mock_dir / f"{accession}.json").open("w") as handle:
        json.dump([record], handle)


def run_downloader(tmpdir, accession, layout, mock_dir, extra=None):
    outdir = tmpdir / "out"
    argv = [
        "--accession",
        accession,
        "--layout",
        layout,
        "--outdir",
        str(outdir),
        "--mock-response-dir",
        str(mock_dir),
        "--max-attempts",
        "2",
        "--retry-wait-seconds",
        "0",
        "--timeout-seconds",
        "1",
    ]
    if extra:
        argv.extend(extra)
    stderr = io.StringIO()
    with redirect_stderr(stderr):
        code = download_sra_fastq.main(argv)
    return code, outdir


def test_paired_download():
    with tempfile.TemporaryDirectory() as temp:
        tmpdir = Path(temp)
        mock_dir = tmpdir / "mock"
        source_dir = tmpdir / "source"
        mock_dir.mkdir()
        source_dir.mkdir()
        r1 = source_dir / "r1.fastq.gz"
        r2 = source_dir / "r2.fastq.gz"
        write_fastq_gz(r1, "SRR000001/1")
        write_fastq_gz(r2, "SRR000001/2")
        write_ena_json(
            mock_dir,
            "SRR000001",
            "PAIRED",
            [r1.as_uri(), r2.as_uri()],
            [md5sum(r1), md5sum(r2)],
            [r1.stat().st_size, r2.stat().st_size],
        )

        code, outdir = run_downloader(tmpdir, "SRR000001", "paired", mock_dir)
        assert code == 0
        assert (outdir / "SRR000001_1.fastq.gz").read_bytes() == r1.read_bytes()
        assert (outdir / "SRR000001_2.fastq.gz").read_bytes() == r2.read_bytes()
        report = json.loads((outdir / "SRR000001.ena_download_report.json").read_text())
        assert report["ena"]["ena_layout"] == "PAIRED"
        assert len(report["downloads"]) == 2


def test_md5_mismatch_fails():
    with tempfile.TemporaryDirectory() as temp:
        tmpdir = Path(temp)
        mock_dir = tmpdir / "mock"
        source_dir = tmpdir / "source"
        mock_dir.mkdir()
        source_dir.mkdir()
        fastq = source_dir / "single.fastq.gz"
        write_fastq_gz(fastq, "SRR000002")
        write_ena_json(mock_dir, "SRR000002", "SINGLE", [fastq.as_uri()], ["0" * 32], [fastq.stat().st_size])

        code, outdir = run_downloader(tmpdir, "SRR000002", "single", mock_dir)
        assert code == 1
        assert not (outdir / "SRR000002.fastq.gz").exists()
        assert not any(outdir.glob("*.part"))
        report = json.loads((outdir / "SRR000002.ena_download_report.json").read_text())
        assert "md5" in report["error"]


def test_missing_source_retries():
    with tempfile.TemporaryDirectory() as temp:
        tmpdir = Path(temp)
        mock_dir = tmpdir / "mock"
        mock_dir.mkdir()
        missing = tmpdir / "missing.fastq.gz"
        write_ena_json(mock_dir, "SRR000003", "SINGLE", [missing.as_uri()])

        code, outdir = run_downloader(tmpdir, "SRR000003", "single", mock_dir, ["--no-verify-md5"])
        assert code == 1
        report = json.loads((outdir / "SRR000003.ena_download_report.json").read_text())
        assert "after 2 attempt" in report["error"]


def test_layout_mismatch_fails():
    with tempfile.TemporaryDirectory() as temp:
        tmpdir = Path(temp)
        mock_dir = tmpdir / "mock"
        source_dir = tmpdir / "source"
        mock_dir.mkdir()
        source_dir.mkdir()
        fastq = source_dir / "single.fastq.gz"
        write_fastq_gz(fastq, "SRR000004")
        write_ena_json(mock_dir, "SRR000004", "SINGLE", [fastq.as_uri()])

        code, outdir = run_downloader(tmpdir, "SRR000004", "paired", mock_dir)
        assert code == 1
        report = json.loads((outdir / "SRR000004.ena_download_report.json").read_text())
        assert "incompatible" in report["error"]


def main():
    failures = 0
    for test in (test_paired_download, test_md5_mismatch_fails, test_missing_source_retries, test_layout_mismatch_fails):
        try:
            test()
        except Exception as exc:
            failures += 1
            print(f"ERROR: {test.__name__}: {exc}", file=sys.stderr)
    if failures:
        return 1
    print("TITAN SRA/ENA download tests OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
