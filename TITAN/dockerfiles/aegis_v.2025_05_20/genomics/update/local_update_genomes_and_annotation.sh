#!/bin/bash

#SBATCH --job-name=update
#SBATCH --output=../logs/updates/update_%j.log
#SBATCH --qos=medium
#SBATCH --time=7-00:00:0
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=4
#SBATCH --mem=50gb

mkdir -p ../logs/updates/

# --ignore-existing allows a faster backup
rsync -avz --ignore-existing --exclude={"*.py","*.pptx","*.pdf"} -e "ssh -p 4144" --progress "tomslab3@147.156.207.102:/media/tomslab3/hard_drive/Documents/10.postdoc_i2sysbio/genomes_and_annotation/" "../../../genomes_and_annotation/"

# --ignore-existing allows a faster backup
rsync -avz --ignore-existing --exclude={"*.py","*.pptx","*.pdf"} -e "ssh -p 4144" --progress "tomslab2@147.156.207.102:/media/tomslab3/hard_drive/Documents/10.postdoc_i2sysbio/genomes_and_annotation/" "../../../genomes_and_annotation/"