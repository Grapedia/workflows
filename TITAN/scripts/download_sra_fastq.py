#!/usr/bin/env python3
import argparse
import hashlib
import json
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


ENA_API_URL = "https://www.ebi.ac.uk/ena/portal/api/filereport"
ENA_FIELDS = ",".join(
    [
        "run_accession",
        "library_layout",
        "fastq_ftp",
        "fastq_md5",
        "fastq_bytes",
        "instrument_platform",
        "library_strategy",
    ]
)
USER_AGENT = "TITAN/1.0 ENA downloader"
TRANSIENT_ERRORS = (
    TimeoutError,
    ConnectionError,
    urllib.error.HTTPError,
    urllib.error.URLError,
)


def ena_report_url(accession):
    query = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": ENA_FIELDS,
            "format": "json",
            "download": "false",
        }
    )
    return f"{ENA_API_URL}?{query}"


def split_field(value):
    return [item for item in (value or "").split(";") if item]


def normalize_url(value):
    if not value:
        return value
    parsed = urllib.parse.urlparse(value)
    if parsed.scheme:
        return value
    if value.startswith("ftp.sra.ebi.ac.uk/") or value.startswith("era-fasp@"):
        return f"ftp://{value}"
    return value


def fetch_ena_records(accession, timeout_seconds, mock_response_dir=None):
    if mock_response_dir:
        path = Path(mock_response_dir) / f"{accession}.json"
        with path.open() as handle:
            return json.load(handle)

    request = urllib.request.Request(
        ena_report_url(accession),
        headers={"Accept": "application/json", "User-Agent": USER_AGENT},
    )
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
        return json.loads(response.read().decode("utf-8"))


def select_record(accession, records):
    if not isinstance(records, list):
        raise ValueError("ENA response is not a JSON list")
    if not records:
        raise ValueError(f"ENA returned no read_run record for {accession}")
    matching = [record for record in records if record.get("run_accession") == accession]
    record = matching[0] if matching else records[0]
    layout = (record.get("library_layout") or "").upper()
    if layout not in {"SINGLE", "PAIRED"}:
        raise ValueError(f"ENA library_layout for {accession} is unsupported: {layout or '<empty>'}")

    urls = [normalize_url(url) for url in split_field(record.get("fastq_ftp"))]
    md5s = split_field(record.get("fastq_md5"))
    sizes = split_field(record.get("fastq_bytes"))
    if not urls:
        raise ValueError(f"ENA record for {accession} does not expose FASTQ URLs")
    if layout == "PAIRED" and len(urls) < 2:
        raise ValueError(f"ENA paired record for {accession} does not expose two FASTQ URLs")
    return {
        "accession": accession,
        "ena_layout": layout,
        "urls": urls,
        "md5s": md5s,
        "sizes": sizes,
        "instrument_platform": record.get("instrument_platform") or "",
        "library_strategy": record.get("library_strategy") or "",
    }


def requested_files(record, requested_layout):
    accession = record["accession"]
    ena_layout = record["ena_layout"]
    if requested_layout == "paired":
        if ena_layout != "PAIRED":
            raise ValueError(f"{accession}: samplesheet layout paired is incompatible with ENA layout {ena_layout}")
        return [
            (record["urls"][0], f"{accession}_1.fastq.gz", optional(record["md5s"], 0), optional(record["sizes"], 0)),
            (record["urls"][1], f"{accession}_2.fastq.gz", optional(record["md5s"], 1), optional(record["sizes"], 1)),
        ]
    if requested_layout in {"single", "long"}:
        if ena_layout != "SINGLE":
            raise ValueError(f"{accession}: samplesheet layout {requested_layout} is incompatible with ENA layout {ena_layout}")
        return [(record["urls"][0], f"{accession}.fastq.gz", optional(record["md5s"], 0), optional(record["sizes"], 0))]
    raise ValueError(f"{accession}: unsupported requested layout {requested_layout}")


def optional(values, index):
    return values[index] if index < len(values) else ""


def md5sum(path):
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def expected_size(value):
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def verify_download(path, expected_md5, expected_bytes, verify_md5):
    if not path.exists() or path.stat().st_size == 0:
        raise ValueError(f"{path.name}: downloaded file is empty or missing")
    size = expected_size(expected_bytes)
    if size is not None and path.stat().st_size != size:
        raise ValueError(f"{path.name}: expected {size} bytes, observed {path.stat().st_size}")
    if verify_md5 and expected_md5:
        observed = md5sum(path)
        if observed.lower() != expected_md5.lower():
            raise ValueError(f"{path.name}: expected md5 {expected_md5}, observed {observed}")


def copy_or_download(url, destination, timeout_seconds):
    parsed = urllib.parse.urlparse(url)
    if parsed.scheme == "file":
        source = Path(urllib.request.url2pathname(parsed.path))
        shutil.copyfile(source, destination)
        return
    if not parsed.scheme:
        shutil.copyfile(Path(url), destination)
        return
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response, destination.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def download_with_retry(url, destination, expected_md5, expected_bytes, args):
    attempts = []
    part = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(1, args.max_attempts + 1):
        if part.exists():
            part.unlink()
        try:
            copy_or_download(url, part, args.timeout_seconds)
            verify_download(part, expected_md5, expected_bytes, args.verify_md5)
            part.replace(destination)
            attempts.append({"attempt": attempt, "status": "success"})
            return attempts
        except TRANSIENT_ERRORS as exc:
            attempts.append({"attempt": attempt, "status": "retryable_error", "error": str(exc)})
        except Exception as exc:
            attempts.append({"attempt": attempt, "status": "error", "error": str(exc)})
        if attempt < args.max_attempts:
            time.sleep(args.retry_wait_seconds)
    if part.exists():
        part.unlink()
    raise RuntimeError(f"failed to download {url} after {args.max_attempts} attempt(s): {attempts[-1]['error']}")


def write_report(path, report):
    with path.open("w") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write("\n")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Download SRA/ENA FASTQ files using ENA API metadata.")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--layout", required=True, choices=["single", "paired", "long"])
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--timeout-seconds", type=int, default=300)
    parser.add_argument("--max-attempts", type=int, default=3)
    parser.add_argument("--retry-wait-seconds", type=float, default=30)
    parser.add_argument("--verify-md5", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--mock-response-dir", type=Path, default=None, help=argparse.SUPPRESS)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    report = {"accession": args.accession, "requested_layout": args.layout, "downloads": []}
    report_path = args.outdir / f"{args.accession}.ena_download_report.json"
    try:
        records = fetch_ena_records(args.accession, args.timeout_seconds, args.mock_response_dir)
        record = select_record(args.accession, records)
        report["ena"] = {key: record[key] for key in ("ena_layout", "instrument_platform", "library_strategy")}
        for url, filename, expected_md5, expected_bytes in requested_files(record, args.layout):
            destination = args.outdir / filename
            attempts = download_with_retry(url, destination, expected_md5, expected_bytes, args)
            report["downloads"].append(
                {
                    "url": url,
                    "path": str(destination),
                    "md5": expected_md5,
                    "bytes": expected_bytes,
                    "attempts": attempts,
                }
            )
        write_report(report_path, report)
        return 0
    except Exception as exc:
        report["error"] = str(exc)
        write_report(report_path, report)
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
