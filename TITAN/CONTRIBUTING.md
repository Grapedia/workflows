# Contributing To TITAN

## Table Of Contents

- [Development References](#development-references)
- [Developer Quality Contract](#developer-quality-contract)
- [Validation And CI](#validation-and-ci)


TITAN uses Nextflow DSL2 modules and subworkflows with a stable public graph.
Development changes should preserve the contracts below.

## Development References

Read these references before changing workflow wiring, module contracts,
container policy or resource labels:

* [Architecture](docs/development/ARCHITECTURE.md)
* [Nextflow DSL2 conventions](docs/development/nextflow-dsl2-conventions.md)
* [Container locks](docs/development/container-locks.md)

## Developer Quality Contract

* Modules should keep workflow logic local to `input:`, `output:`, `script:`
  and `stub:` blocks. Domain options should be passed as `val` or `path`
  inputs instead of being read directly from `params`; `params` remains
  acceptable for containers, labels and publish locations.
* Existing process names are kept stable to avoid breaking `withName`
  selectors, traces and resume behavior. New processes should use lower snake
  case matching the module filename.
* Every process must use an explicit label from the resource policy in
  [conf/base.config](conf/base.config).
* Process outputs should be named with `emit:` and should use `path(...)` for
  files. Avoid broad output globs that can capture temporary files.
* Every process should emit `versions.yml`. For tools where a reliable runtime
  `--version` command is not available, record at least the process identity
  and configured container/tool version.
* Shell scripts in modules must use `set -euo pipefail`. Complex repeated
  shell workflows should live under `scripts/` and be passed into modules as
  `path` inputs.
* Workflow closures must not read task output files with `file(...).text`;
  emit values from the producing process or move parsing into a process/script.
* User-facing parameters should be documented and typed in
  [nextflow_schema.json](nextflow_schema.json).
* Keep `scripts/run-tests.sh` passing. Add focused tests or static checks when
  changing tuple contracts, module outputs, labels, publication behavior or
  shared helper scripts.

`nf-test` specs are the preferred direction for future module-level tests once
contracts are stable. Until then, TITAN uses Python/static checks plus the full
`test` profile stub run in `scripts/run-tests.sh`.

## Validation And CI

Local quick validation:

```bash
scripts/run-tests.sh
```

Targeted validation:

```bash
python3 scripts/validate_container_pins.py
python3 scripts/validate_profiles.py
python3 scripts/validate_minimal_test_data.py
python3 scripts/test_validate_inputs.py
python3 scripts/test_validate_final_annotation.py
nextflow run main.nf -profile test -stub-run -ansi-log false
```

GitHub Actions runs the same quick suite when files under `TITAN/**` or the
TITAN CI workflow change.
