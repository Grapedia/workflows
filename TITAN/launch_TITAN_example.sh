#!/usr/bin/env bash
# Exit immediately if a command exits with a non-zero status
# Ensures AEGIS doesn't run if generate_evidence_data fails
set -e
# Navigate to the project workflow directory
cd /home/avelt/data2/2024_assembly_GW_RI_hifiasm/Riesling/2025_genes_annotation/workflows/TITAN
# Load required Nextflow module
module load nextflow/24.04.3
# Run the 'generate_evidence_data' workflow and generate its DAG
# nextflow run main.nf \
#  -with-dag dag_evidence_data.png \
#  --workflow generate_evidence_data`
# Run the 'aegis' workflow and generate its DAG
nextflow run main.nf \
  -with-dag dag_aegis.png \
  --workflow aegis
