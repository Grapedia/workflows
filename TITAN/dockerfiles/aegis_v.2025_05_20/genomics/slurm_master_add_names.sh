#!/bin/bash

#SBATCH --job-name=genomics
#SBATCH -o logs/annotation_objects/annotation_objects_manager_%j_out.log
#SBATCH -e logs/annotation_objects/annotation_objects_manager_%j_err.log
#SBATCH --ntasks=1    
#SBATCH --time=0-10:00:0       
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=2
#SBATCH --qos=short

module load python/3.11.4

python3 -u annotation_v5.1_add_gene_names.py