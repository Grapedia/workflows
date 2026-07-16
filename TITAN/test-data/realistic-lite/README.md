# TITAN Realistic Lite Fixture

This fixture is a tiny synthetic end-to-end dataset for non-stub TITAN smoke runs.
It contains a 3.6 kb target assembly, a nearly identical reference assembly, one
GFF3 gene model, local single/paired/long RNA-seq reads and one protein evidence
FASTA.

Run with:

```bash
nextflow -c test-data/realistic-lite/real_local.config run main.nf -profile local -w test-work/realistic-lite -ansi-log false
```

EGAPx still requires a prepared local cache. The config points `egapx_runner_dir`
to `${projectDir}/.egapx_runner` and `egapx_local_cache_dir` to
`${projectDir}/.egapx_cache`.
