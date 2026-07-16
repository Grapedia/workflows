#!/usr/bin/env python3
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


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


def main():
    result = validate_no_container_options()
    if result:
        return result

    expectations = {
        "test": ["executor = 'local'", "enabled = false", "workDir ="],
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

    print("TITAN profiles OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
