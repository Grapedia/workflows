#!/bin/bash

# This script launches DIAMOND to create databases and run BLASTp searches.
# Usage : ./diamond_script.sh -q <proteins_fasta_query_file> -t <threads> -d <proteins_fasta_list_reference> -o <output_dir>

# default value
THREADS=5

# Usage function
print_usage() {
    echo "Usage : $0 -q <query_file> -t <threads> -d <databases> -o <output_dir>"
    echo "    -q    Path to the FASTA file of unique proteins for the assembly to annotate"
    echo "    -t    Number of threads (default : 10)"
    echo "    -o    Output directory"
    echo "    -d    Absolute path to proteins fasta files to be used as databases, separated by commas (example: /path/to/arabidopsis.proteins.fasta,/path/to/viridiplantae.proteins.fasta,/path/to/eudicotyledons.proteins.fasta)"
    exit 1
}

# Reading the arguments
while getopts "q:t:d:o:h" option; do
    case "$option" in
        q) QUERY="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        d) DATABASES="$OPTARG" ;;
		o) OUTDIR="$OPTARG" ;;
        h) print_usage ;;
        *) print_usage ;;
    esac
done

# Checking arguments
if [ -z "$QUERY" ] || [ -z "$DATABASES" ] || [ -z "$OUTDIR" ]; then
    echo "Error: the query file, output directory and databases must be specified with -q, -o and -d"
    print_usage
fi

echo "------------------------------ QUERY FILE LOADED :"
echo "QUERY: $QUERY"

# IFS=‘,’ database loop
IFS=',' read -ra DB_ARRAY <<< "$DATABASES"
for DB in "${DB_ARRAY[@]}"; do
    DB_PATH="${DB}"
    BASENAME=$( basename "$DB_PATH" )
    OUTPUT_PATH="$OUTDIR/${BASENAME}_vs_assembly.diamond"

    echo "Processing database: $DB_PATH"

    # Creating the DIAMOND database
    diamond makedb --threads "$THREADS" --in "$DB_PATH" \
    -d "$OUTDIR/${BASENAME}"

    # Run DIAMOND BLASTp
    diamond blastp --threads "$THREADS" --db "$OUTDIR/${BASENAME}" \
    --ultra-sensitive --out "$OUTPUT_PATH" \
    --outfmt 6 --query "$QUERY" --max-target-seqs 1 --evalue 1e-3
done
