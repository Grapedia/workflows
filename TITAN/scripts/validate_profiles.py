#!/usr/bin/env python3
import subprocess
import sys
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESOURCE_LABELS = {
    "process_low",
    "process_trim",
    "process_index",
    "process_alignment",
    "process_transcriptome",
    "process_prediction",
    "process_merge",
    "process_aegis",
    "process_rfam",
}


def fail(message):
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def run_config(profile):
    return subprocess.check_output(
        ["nextflow", "config", "-profile", profile],
        cwd=ROOT,
        text=True,
        stderr=subprocess.STDOUT,
    )


def require_contains(profile, config, expected):
    missing = [item for item in expected if item not in config]
    if missing:
        return fail(f"profile {profile} is missing expected config fragment(s): {', '.join(missing)}")
    return 0


def validate_no_container_options():
    offenders = []
    for path in sorted((ROOT / "modules").glob("*.nf")):
        for line_number, line in enumerate(path.read_text().splitlines(), 1):
            if "containerOptions" in line:
                offenders.append(f"{path.relative_to(ROOT)}:{line_number}: {line.strip()}")
    if offenders:
        return fail("module containerOptions are not portable to Apptainer:\n" + "\n".join(offenders))
    return 0


def validate_resource_labels():
    base_config = (ROOT / "conf" / "base.config").read_text()
    missing_config = [label for label in sorted(RESOURCE_LABELS) if f"withLabel: {label}" not in base_config]
    if missing_config:
        return fail("conf/base.config is missing resource label(s): " + ", ".join(missing_config))

    offenders = []
    missing_labels = []
    unknown_labels = []
    for path in sorted((ROOT / "modules").glob("*.nf")):
        text = path.read_text()
        labels = re.findall(r"""label\s+['"](process_[^'"]+)['"]""", text)
        if not labels:
            missing_labels.append(str(path.relative_to(ROOT)))
        unknown_labels.extend(
            f"{path.relative_to(ROOT)}: {label}"
            for label in labels
            if label not in RESOURCE_LABELS
        )
        for line_number, line in enumerate(text.splitlines(), 1):
            stripped = line.strip()
            if stripped.startswith("cpus "):
                offenders.append(f"{path.relative_to(ROOT)}:{line_number}: {stripped}")
    if missing_labels:
        return fail("module process resource labels are missing:\n" + "\n".join(missing_labels))
    if unknown_labels:
        return fail("unknown module process resource label(s):\n" + "\n".join(unknown_labels))
    if offenders:
        return fail("module cpus directives must be configured through labels or withName selectors:\n" + "\n".join(offenders))
    return 0


def validate_configured_cpus(profile, config):
    if "cpus = null" in config:
        return fail(f"profile {profile} resolved at least one process CPU directive to null")
    return 0


def main():
    result = validate_no_container_options()
    if result:
        return result
    result = validate_resource_labels()
    if result:
        return result

    expectations = {
        "test": ["executor = 'local'", "enabled = false", "workDir =", "withLabel:process_prediction"],
        "local": ["executor = 'local'", "docker", "enabled = true", "egapx_executor = 'docker'"],
        "apptainer,test": ["apptainer", "enabled = true", "egapx_executor = 'singularity'"],
        "test,slurm,apptainer": ["executor = 'slurm'", "queueSize", "submitRateLimit", "egapx_executor = 'singularity'"],
    }

    for profile, expected in expectations.items():
        try:
            config = run_config(profile)
        except subprocess.CalledProcessError as exc:
            print(exc.output, file=sys.stderr)
            return fail(f"nextflow config failed for profile {profile}")
        result = require_contains(profile, config, expected)
        if result:
            return result
        result = validate_configured_cpus(profile, config)
        if result:
            return result

    print("TITAN profiles OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
