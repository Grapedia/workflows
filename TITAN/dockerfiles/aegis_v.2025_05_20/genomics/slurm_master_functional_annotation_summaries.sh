#!/bin/bash

#SBATCH --job-name=fun_a
#SBATCH -o logs/functional_annotation/fun_a_%j_out.log
#SBATCH -e logs/functional_annotation/fun_a_%j_err.log
#SBATCH --ntasks=1    
#SBATCH --time=2-00:00:00       
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=2
#SBATCH --qos=medium

module load python/3.11.4

python3 -u formatting_functional_annotation_all_species.py

