#!/bin/bash

tail -n +2 "$1" | while IFS=',' read -r col1 col2; do
    echo "$col1 /protein_path/$col2"
done
