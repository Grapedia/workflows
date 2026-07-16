# TITAN test strategy

`scripts/run-tests.sh` is the local quick test entry point. It is intentionally dependency-light and uses only Python, Bash and Nextflow so it can run on a development machine without Slurm, Docker or Apptainer execution.

The quick suite covers:

* static container pin checks;
* profile resolution and portability checks;
* minimal fixture validation;
* unit-style input validator cases;
* full TITAN `test` profile execution in Nextflow stub mode;
* one negative Nextflow case proving invalid inputs fail before heavy execution.

Future `nf-test` specs should live under this directory and target individual modules or subworkflows once their contracts are stable. Keep those tests on the synthetic fixtures under `test-data/minimal`.
