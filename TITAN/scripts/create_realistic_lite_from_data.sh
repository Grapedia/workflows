#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
OUT_DIR="${1:-${ROOT_DIR}/data_example/realistic-lite}"
SOURCE_DATA_DIR="${SOURCE_DATA_DIR:-${ROOT_DIR}/data}"

FASTQ_RECORDS="${FASTQ_RECORDS:-50000}"
LONG_RECORDS="${LONG_RECORDS:-3000}"
PROTEIN_RECORDS="${PROTEIN_RECORDS:-5000}"
CHROMOSOME="${CHROMOSOME:-chr01}"

PREVIOUS_ASSEMBLY="${SOURCE_DATA_DIR}/assemblies/v4_genome_ref.fasta"
NEW_ASSEMBLY="${SOURCE_DATA_DIR}/assemblies/T2T_ref.fasta"
PREVIOUS_ANNOTATIONS="${SOURCE_DATA_DIR}/annotations/v4_3_just_ref.gff3"
PAIRED_SAMPLE="depleted-stranded_Vitis14_S16"
SINGLE_SAMPLE="PolyA_Vitis02_S1"
LONG_SAMPLE="PN40024_IsoSeq_chr01"
PAIRED_R1="${SOURCE_DATA_DIR}/RNAseq_data/${PAIRED_SAMPLE}_1.fastq.gz"
PAIRED_R2="${SOURCE_DATA_DIR}/RNAseq_data/${PAIRED_SAMPLE}_2.fastq.gz"
SINGLE_R1="${SOURCE_DATA_DIR}/RNAseq_data/${SINGLE_SAMPLE}_1.fastq.gz"
LONG_FASTA="${SOURCE_DATA_DIR}/RNAseq_data/hq_transcripts.RI_rmv_Antonio.fasta"
LONG_GMAP="${SOURCE_DATA_DIR}/RNAseq_data/IsoSeq/gmap.spliced_alignments_ChangedChrName.gff3"
PROTEIN_SWISSPROT="${SOURCE_DATA_DIR}/protein_data/Viridiplantae_swissprot.cleaned.fasta"
PROTEIN_ORTHODB="${SOURCE_DATA_DIR}/protein_data/eudicotyledons_odb10.cleaned.fasta"

require_file() {
  local path="$1"
  if [[ ! -s "$path" ]]; then
    echo "ERROR: required source file is missing or empty: $path" >&2
    exit 1
  fi
}

gzip_cmd() {
  if command -v pigz >/dev/null 2>&1; then
    pigz -c
  else
    gzip -c
  fi
}

extract_fasta_sequence() {
  local seq_id="$1"
  local input="$2"
  local output="$3"
  awk -v seq_id="$seq_id" '
    /^>/ {
      current = $1
      sub(/^>/, "", current)
      keep = current == seq_id
    }
    keep { print }
  ' "$input" > "$output"
  [[ -s "$output" ]] || {
    echo "ERROR: sequence ${seq_id} not found in ${input}" >&2
    exit 1
  }
}

extract_first_fasta_records() {
  local max_records="$1"
  local input="$2"
  local output="$3"
  awk -v max_records="$max_records" '
    /^>/ {
      record_count++
      keep = record_count <= max_records
    }
    keep { print }
  ' "$input" > "$output"
  [[ -s "$output" ]] || {
    echo "ERROR: no FASTA records written to ${output}" >&2
    exit 1
  }
}

extract_fasta_records_by_id() {
  local ids="$1"
  local input="$2"
  local output="$3"
  awk '
    FNR == NR {
      wanted[$1] = 1
      next
    }
    /^>/ {
      current = $1
      sub(/^>/, "", current)
      keep = current in wanted
    }
    keep { print }
  ' "$ids" "$input" > "$output"
  [[ -s "$output" ]] || {
    echo "ERROR: no selected FASTA records written to ${output}" >&2
    exit 1
  }
}

subset_fastq_gz() {
  local records="$1"
  local input="$2"
  local output="$3"
  zcat "$input" | awk -v records="$records" 'NR <= records * 4 { print }' | gzip_cmd > "$output"
  [[ -s "$output" ]] || {
    echo "ERROR: no FASTQ records written to ${output}" >&2
    exit 1
  }
}

for source in \
  "$PREVIOUS_ASSEMBLY" \
  "$NEW_ASSEMBLY" \
  "$PREVIOUS_ANNOTATIONS" \
  "$PAIRED_R1" \
  "$PAIRED_R2" \
  "$SINGLE_R1" \
  "$LONG_FASTA" \
  "$LONG_GMAP" \
  "$PROTEIN_SWISSPROT" \
  "$PROTEIN_ORTHODB"
do
  require_file "$source"
done

mkdir -p \
  "${OUT_DIR}/assemblies" \
  "${OUT_DIR}/annotations" \
  "${OUT_DIR}/rnaseq" \
  "${OUT_DIR}/proteins" \
  "${OUT_DIR}/egapx" \
  "${OUT_DIR}/egapx_config" \
  "${OUT_DIR}/logs"

echo "Creating realistic lite fixture in ${OUT_DIR}"
echo "Chromosome: ${CHROMOSOME}"
echo "FASTQ records per short-read file: ${FASTQ_RECORDS}"
echo "Long-read FASTA records: ${LONG_RECORDS}"
echo "Protein records per FASTA: ${PROTEIN_RECORDS}"

extract_fasta_sequence "$CHROMOSOME" "$PREVIOUS_ASSEMBLY" "${OUT_DIR}/assemblies/v4_genome_ref.${CHROMOSOME}.fasta"
extract_fasta_sequence "$CHROMOSOME" "$NEW_ASSEMBLY" "${OUT_DIR}/assemblies/T2T_ref.${CHROMOSOME}.fasta"

awk -v chrom="$CHROMOSOME" '
  BEGIN { print "##gff-version 3" }
  !/^#/ && $1 == chrom { print }
' "$PREVIOUS_ANNOTATIONS" > "${OUT_DIR}/annotations/v4_3_just_ref.${CHROMOSOME}.gff3"

subset_fastq_gz "$FASTQ_RECORDS" "$PAIRED_R1" "${OUT_DIR}/rnaseq/${PAIRED_SAMPLE}_1.fastq.gz"
subset_fastq_gz "$FASTQ_RECORDS" "$PAIRED_R2" "${OUT_DIR}/rnaseq/${PAIRED_SAMPLE}_2.fastq.gz"
subset_fastq_gz "$FASTQ_RECORDS" "$SINGLE_R1" "${OUT_DIR}/rnaseq/${SINGLE_SAMPLE}.fastq.gz"

awk -v chrom="$CHROMOSOME" -v max_records="$LONG_RECORDS" '
  BEGIN { FS = "\t" }
  $1 == chrom && $3 == "cDNA_match" {
    split($9, attrs, ";")
    for (i in attrs) {
      if (attrs[i] ~ /^Target=/) {
        target = attrs[i]
        sub(/^Target=/, "", target)
        split(target, parts, " ")
        id_number = parts[1]
        if (id_number ~ /^Ware_HQ_transcript_[0-9]+_/) {
          sub(/^Ware_HQ_transcript_/, "", id_number)
          sub(/_.*/, "", id_number)
          id = "Ware_HQ_transcript/" id_number
          if (!seen[id]++) {
            print id
            selected++
            if (selected >= max_records) {
              exit
            }
          }
        }
      }
    }
  }
' "$LONG_GMAP" > "${OUT_DIR}/logs/${LONG_SAMPLE}.selected_ids.txt"

extract_fasta_records_by_id \
  "${OUT_DIR}/logs/${LONG_SAMPLE}.selected_ids.txt" \
  "$LONG_FASTA" \
  "${OUT_DIR}/rnaseq/${LONG_SAMPLE}.fasta"

extract_first_fasta_records "$PROTEIN_RECORDS" "$PROTEIN_SWISSPROT" "${OUT_DIR}/proteins/Viridiplantae_swissprot.cleaned.${PROTEIN_RECORDS}.fasta"
extract_first_fasta_records "$PROTEIN_RECORDS" "$PROTEIN_ORTHODB" "${OUT_DIR}/proteins/eudicotyledons_odb10.cleaned.${PROTEIN_RECORDS}.fasta"

cat > "${OUT_DIR}/RNAseq_samplesheet.csv" <<EOF
sampleID,SRA_or_FASTQ,library_layout
${PAIRED_SAMPLE},FASTQ,paired
${SINGLE_SAMPLE},FASTQ,single
${LONG_SAMPLE},FASTA,long
EOF

cat > "${OUT_DIR}/protein_samplesheet.csv" <<EOF
organism,filename
viridiplantae,data_example/realistic-lite/proteins/Viridiplantae_swissprot.cleaned.${PROTEIN_RECORDS}.fasta
eudicotyledones_orthoDB,data_example/realistic-lite/proteins/eudicotyledons_odb10.cleaned.${PROTEIN_RECORDS}.fasta
EOF

cat > "${OUT_DIR}/egapx/input_egapx.yaml" <<EOF
genome: ${OUT_DIR}/assemblies/T2T_ref.${CHROMOSOME}.fasta
taxid: 29760
organism: Vitis vinifera
annotation_provider: egapx_ncbi
annotation_name_prefix: PN40024_${CHROMOSOME}_lite
locus_tag_prefix: VitviLite
short_reads:
  - - ${PAIRED_SAMPLE}
    - - ${OUT_DIR}/rnaseq/${PAIRED_SAMPLE}_1.fastq.gz
      - ${OUT_DIR}/rnaseq/${PAIRED_SAMPLE}_2.fastq.gz
  - - ${SINGLE_SAMPLE}
    - - ${OUT_DIR}/rnaseq/${SINGLE_SAMPLE}.fastq.gz
long_reads:
  - - ${LONG_SAMPLE}
    - - ${OUT_DIR}/rnaseq/${LONG_SAMPLE}.fasta
EOF

cat > "${OUT_DIR}/egapx_config/process_resources.config" <<'EOF'
params.threads = 7
params.nodes = 1
params.num_cpus_per_node = 8

process {
    memory = 56.GB
    cpus = 7
    ext.align_mem = '13'
    time = 24.h
    errorStrategy = 'retry'
    maxRetries = 1

    withLabel: 'long_job' {
        time = 36.h
    }
    withLabel: 'small_mem' {
        memory = 8.GB
    }
    withLabel: 'med_mem' {
        memory = 32.GB
    }
    withLabel: 'large_mem' {
        memory = 56.GB
    }
    withLabel: 'single_cpu' {
        ext {
            split_jobs = 1
            threads = 1
        }
        cpus = 1
    }
    withLabel: 'multi_cpu' {
        ext {
            split_jobs = 1
            threads = params.threads
        }
        cpus = params.threads
    }
    withLabel: 'multi_node' {
        ext {
            split_jobs = params.nodes
            threads = params.threads
        }
        cpus = params.threads
    }
    withLabel: 'gpx_submitter' {
        ext {
            split_jobs = params.nodes
            threads = 1
        }
        cpus = 1
    }
}
EOF

cat > "${OUT_DIR}/egapx_config/singularity.config" <<'EOF'
includeConfig "./process_resources.config"

singularity {
    enabled = true
    autoMounts = true
    cacheDir = System.getenv('SINGULARITY_CACHEDIR') ?: System.getenv('APPTAINER_CACHEDIR') ?: './singularity-cache'
    envWhitelist = 'APPTAINER_CACHEDIR,SINGULARITY_CACHEDIR,APPTAINER_TMPDIR,SINGULARITY_TMPDIR,PYTHONNOUSERSITE'
}

env {
    APPTAINER_CACHEDIR = System.getenv('APPTAINER_CACHEDIR') ?: System.getenv('SINGULARITY_CACHEDIR') ?: './singularity-cache'
    SINGULARITY_CACHEDIR = System.getenv('SINGULARITY_CACHEDIR') ?: System.getenv('APPTAINER_CACHEDIR') ?: './singularity-cache'
    APPTAINER_TMPDIR = System.getenv('APPTAINER_TMPDIR') ?: System.getenv('SINGULARITY_TMPDIR') ?: './apptainer-tmp'
    SINGULARITY_TMPDIR = System.getenv('SINGULARITY_TMPDIR') ?: System.getenv('APPTAINER_TMPDIR') ?: './apptainer-tmp'
    PYTHONNOUSERSITE = '1'
    DEBUGINFOD_URLS = '/dev/null'
}

process {
    executor = 'local'
    cache = 'lenient'
}
EOF

cat > "${OUT_DIR}/real_slurm_apptainer.config" <<EOF
env {
  APPTAINER_CACHEDIR = "\${projectDir}/.apptainer-cache"
  SINGULARITY_CACHEDIR = "\${projectDir}/.apptainer-cache"
  APPTAINER_TMPDIR = "\${projectDir}/.apptainer-tmp"
  SINGULARITY_TMPDIR = "\${projectDir}/.apptainer-tmp"
}

apptainer {
  envWhitelist = 'APPTAINER_CACHEDIR,SINGULARITY_CACHEDIR,APPTAINER_TMPDIR,SINGULARITY_TMPDIR'
}

singularity {
  envWhitelist = 'APPTAINER_CACHEDIR,SINGULARITY_CACHEDIR,APPTAINER_TMPDIR,SINGULARITY_TMPDIR'
}

params {
  output_dir = "\${projectDir}/test-results/data-example-realistic-lite"
  previous_assembly = "\${projectDir}/data_example/realistic-lite/assemblies/v4_genome_ref.${CHROMOSOME}.fasta"
  new_assembly = "\${projectDir}/data_example/realistic-lite/assemblies/T2T_ref.${CHROMOSOME}.fasta"
  previous_annotations = "\${projectDir}/data_example/realistic-lite/annotations/v4_3_just_ref.${CHROMOSOME}.gff3"
  RNAseq_samplesheet = "\${projectDir}/data_example/realistic-lite/RNAseq_samplesheet.csv"
  RNAseq_data_dir = "\${projectDir}/data_example/realistic-lite/rnaseq"
  protein_samplesheet = "\${projectDir}/data_example/realistic-lite/protein_samplesheet.csv"
  egapx_paramfile = "\${projectDir}/data_example/realistic-lite/egapx/input_egapx.yaml"
  edta_cpus = 8
  egapx_cpus = 8
  egapx_python = "/trinity/shared/apps/python_3.12/bin/python3"
  egapx_data_version = "egapxsupportdata_20251017"
  egapx_config_dir = "\${projectDir}/data_example/realistic-lite/egapx_config"
  diamond2go_cpus = 4
  ena_retry_wait_seconds = 0
  publish_intermediates = true
  slurm_export_env = "--export=ALL,APPTAINER_CACHEDIR=\${projectDir}/.apptainer-cache,SINGULARITY_CACHEDIR=\${projectDir}/.apptainer-cache,APPTAINER_TMPDIR=\${projectDir}/.apptainer-tmp,SINGULARITY_TMPDIR=\${projectDir}/.apptainer-tmp"
}

process {
  clusterOptions = "--export=ALL,APPTAINER_CACHEDIR=\${projectDir}/.apptainer-cache,SINGULARITY_CACHEDIR=\${projectDir}/.apptainer-cache,APPTAINER_TMPDIR=\${projectDir}/.apptainer-tmp,SINGULARITY_TMPDIR=\${projectDir}/.apptainer-tmp"
  withLabel: process_low { cpus = 1; memory = 4.GB; time = 2.h }
  withLabel: process_medium { cpus = 4; memory = 16.GB; time = 6.h }
  withLabel: process_high { cpus = 8; memory = 48.GB; time = 24.h }
  withLabel: process_index { cpus = 8; memory = 32.GB; time = 8.h }
  withLabel: process_alignment { cpus = 8; memory = 32.GB; time = 8.h }
  withLabel: process_transcriptome { cpus = 8; memory = 32.GB; time = 8.h }
  withLabel: process_prediction { cpus = 8; memory = 64.GB; time = 48.h }
  withLabel: process_merge { cpus = 4; memory = 16.GB; time = 8.h }
  withLabel: process_aegis { cpus = 4; memory = 32.GB; time = 12.h }
}
EOF

cat > "${OUT_DIR}/README.md" <<EOF
# TITAN realistic lite fixture from PN40024 data

This fixture is generated from the real \`data/\` inputs, reduced to \`${CHROMOSOME}\`
and a small but non-trivial RNA/protein evidence set.

Contents:

- target assembly: \`assemblies/T2T_ref.${CHROMOSOME}.fasta\`
- previous assembly: \`assemblies/v4_genome_ref.${CHROMOSOME}.fasta\`
- previous GFF3 annotation: \`annotations/v4_3_just_ref.${CHROMOSOME}.gff3\`
- paired short reads: \`${PAIRED_SAMPLE}_1/2.fastq.gz\`
- single short reads: \`${SINGLE_SAMPLE}.fastq.gz\`
- long reads: \`${LONG_SAMPLE}.fasta\`
- reduced protein evidence: \`${PROTEIN_RECORDS}\` records from each cleaned source FASTA

Regenerate from the repository root with:

\`\`\`bash
FASTQ_RECORDS=${FASTQ_RECORDS} LONG_RECORDS=${LONG_RECORDS} PROTEIN_RECORDS=${PROTEIN_RECORDS} \\
  scripts/create_realistic_lite_from_data.sh
\`\`\`

Validate inputs with:

\`\`\`bash
python3 scripts/validate_inputs.py \\
  --project-dir . \\
  --new-assembly data_example/realistic-lite/assemblies/T2T_ref.${CHROMOSOME}.fasta \\
  --previous-assembly data_example/realistic-lite/assemblies/v4_genome_ref.${CHROMOSOME}.fasta \\
  --previous-annotations data_example/realistic-lite/annotations/v4_3_just_ref.${CHROMOSOME}.gff3 \\
  --rnaseq-samplesheet data_example/realistic-lite/RNAseq_samplesheet.csv \\
  --rnaseq-data-dir data_example/realistic-lite/rnaseq \\
  --protein-samplesheet data_example/realistic-lite/protein_samplesheet.csv \\
  --egapx-paramfile data_example/realistic-lite/egapx/input_egapx.yaml \\
  --egapx-executor apptainer
\`\`\`

Run a cluster smoke test first with \`-stub-run\`, then remove \`-stub-run\` for
the real end-to-end run:

\`\`\`bash
nextflow -c data_example/realistic-lite/real_slurm_apptainer.config \\
  run main.nf -profile slurm,apptainer -stub-run -ansi-log false
\`\`\`
EOF

find "$OUT_DIR" -type f -printf '%s\t%p\n' | sort -n > "${OUT_DIR}/logs/file_manifest.tsv"

echo "Fixture creation complete."
echo "Manifest: ${OUT_DIR}/logs/file_manifest.tsv"
