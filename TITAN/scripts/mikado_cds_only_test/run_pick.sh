#!/usr/bin/env bash
#SBATCH -J mik_cdsonly --nodes=1 -c 16 --mem=64G -t 06:00:00
#SBATCH -o %x-%j.log
# Re-run only `mikado pick` of a production run with pick.clustering.cds_only: true
# (coding transcripts clustered on CDS overlap only, so a gene sitting in a neighbour's UTR is no longer
# merged into its locus). Same database, prepared GTF, genome and configuration as production.
#
#   sbatch --exclude=calcul scripts/mikado_cds_only_test/run_pick.sh <prod_pick_workdir> <out_dir>
#
# <prod_pick_workdir>: Nextflow work dir of the production mikado_pick task (holds mikado.db,
# mikado_configuration.yaml, mikado_prepared.gtf and the genome as symlinks).
# Then compare with scripts/mikado_cds_only_test/compare_loci.py.
set -Eeuo pipefail
SRC=$(readlink -f "${1:?usage: run_pick.sh <prod_pick_workdir> <out_dir>}")
OUT=${2:?usage: run_pick.sh <prod_pick_workdir> <out_dir>}
PROJECT_DIR="${TITAN_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)}}"
IMG="$PROJECT_DIR/.apptainer-cache/quay.io-biocontainers-mikado@sha256-dd6f5a2a2d7fdbab73c835cd0f49bd1444ecaddf8e4cd96fbf0fe24f5ecf5f22.img"

mkdir -p "$OUT"
cd "$OUT"
for f in mikado.db mikado_prepared.gtf; do ln -sf "$(readlink -f "$SRC/$f")" "$f"; done
genome=$(awk '/^reference:/ {r=1} r && /genome:/ {print $2; exit}' "$SRC/mikado_configuration.yaml")
ln -sf "$(readlink -f "$SRC/$genome")" "$genome"
awk '{print} /^pick:$/ {print "  clustering:\n    cds_only: true"}' "$SRC/mikado_configuration.yaml" > mikado_configuration.yaml
grep -q "cds_only: true" mikado_configuration.yaml

export TMPDIR="$PROJECT_DIR/.tmp"
singularity exec --cleanenv --bind /data2 "$IMG" mikado pick --json-conf mikado_configuration.yaml \
  --subloci-out mikado.subloci.gff3 --loci-out mikado.loci.gff3 --procs "${SLURM_CPUS_PER_TASK:-16}"
test -s mikado.loci.gff3 && echo DONE
