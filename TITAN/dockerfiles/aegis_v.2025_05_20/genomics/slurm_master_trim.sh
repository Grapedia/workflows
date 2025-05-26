#!/bin/bash

#SBATCH --job-name=CDS_trim
#SBATCH -o logs/other/CDS_trim_%j_out.log
#SBATCH -e logs/other/CDS_trim_%j_err.log
#SBATCH --ntasks=1    
#SBATCH --time=2-00:00:00       
#SBATCH --mem=20gb
#SBATCH --cpus-per-task=2
#SBATCH --qos=medium

module load python/3.11.4


python3 -u annotation_just_one_object_unofficial_CDS_trim.py

