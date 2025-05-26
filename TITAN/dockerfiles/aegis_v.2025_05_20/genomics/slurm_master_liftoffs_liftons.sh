#!/bin/bash

#SBATCH --job-name=liftoff
#SBATCH -o logs/liftoffs/liftoffs_liftons_%j_out.log
#SBATCH -e logs/liftoffs/liftoffs_liftons_%j_err.log
#SBATCH --ntasks=1    
#SBATCH --time=1-00:00:00       
#SBATCH --mem=30gb
#SBATCH --cpus-per-task=8
#SBATCH --qos=short

module load python/3.11.4
module load singularity/3.11.1

python3 -u PN40024_liftoffs_liftons_snapdragon.py

