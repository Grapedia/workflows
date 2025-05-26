#!/bin/bash

#SBATCH --job-name=genomics   
#SBATCH --output=logs/other/5.1.1_%j.log 
#SBATCH --ntasks=1    
#SBATCH --time=1-00:00:0       
#SBATCH --mem=40gb
#SBATCH --cpus-per-task=4
#SBATCH --qos=medium

module load python/3.11.4

python3 -u annotation_just_one_object_release_5.1.1.py