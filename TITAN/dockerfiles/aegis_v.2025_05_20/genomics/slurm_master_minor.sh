#!/bin/bash

#SBATCH --job-name=genomics
#SBATCH -o logs/other/minor_%j_out.log
#SBATCH -e logs/other/minor_%j_err.log
#SBATCH --ntasks=1    
#SBATCH --time=10:00:0       
#SBATCH --mem=50gb
#SBATCH --cpus-per-task=2
#SBATCH --qos=short

module load python/3.11.4

python3 -u annotation_just_one_object_mark_rRNAs.py