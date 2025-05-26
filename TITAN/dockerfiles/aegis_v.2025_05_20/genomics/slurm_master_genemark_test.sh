#!/bin/bash

#SBATCH --job-name=genomics   
#SBATCH --output=logs/other/genemark_test_%j.log 
#SBATCH --ntasks=1    
#SBATCH --time=1-00:00:0       
#SBATCH --mem=10gb
#SBATCH --cpus-per-task=4
#SBATCH --qos=medium

module load python/3.11.4

python3 -u genemark_gff_test.py