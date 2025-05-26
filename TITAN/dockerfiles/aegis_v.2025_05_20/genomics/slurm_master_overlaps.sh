#!/bin/bash

#SBATCH --job-name=overlaps
#SBATCH -o logs/overlaps/overlaps_%j_out.log
#SBATCH -e logs/overlaps/overlaps_%j_err.log
#SBATCH --ntasks=1    
#SBATCH --time=5-00:00:00       
#SBATCH --mem=140gb
#SBATCH --cpus-per-task=2
#SBATCH --qos=medium

module load python/3.11.4

python3 -u PN40024_overlaps_T2T.py
