# TITAN Nextflow DSL2 conventions

These conventions are intentionally conservative for the current codebase.

## Process naming

Existing process names are kept stable to preserve `withName` selectors, traces and resume behavior. New processes should use lower snake case matching the module filename, for example `eggnog_mapper` in `modules/eggnog_mapper.nf`.

Do not rename existing processes only for style. Rename only when there is a migration plan for selectors, tests and published documentation.

## Module contracts

- Modules may use `params` for containers, labels and publish locations.
- Domain options should be passed as `val` or `path` inputs so task hashes record the effective contract.
- Every module must emit `versions.yml`.
- Shell scripts in modules must start with `set -euo pipefail`.
- Workflow closures must not read staged task files with `file(...).text`; emit values from processes or move parsing into a process/script.

## Shared scripts

When two modules need the same multi-step shell workflow, put the implementation under `scripts/` and pass it as a `path` input. The module should remain responsible for the Nextflow contract and the script should own the command-line workflow.
