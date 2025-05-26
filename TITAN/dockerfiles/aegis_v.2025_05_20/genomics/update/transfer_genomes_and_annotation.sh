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
# --delete allows to delete things in destination, but don't use it if you don't want to lose index genomes or other additional destination files

# echo transferring to Luis

# rsync -avz --delete --exclude={"*.py","*.pptx","*.pdf","*.fai","*.mmi"} -e "ssh -p 4144" --progress "../../../genomes_and_annotation/" "tomslab@147.156.207.134:/media/tomslab/hard_drive/Luis/genomes_and_annotation/"

# echo transfering to raspberry

# rsync -avz --delete --exclude={"*.py","*.pptx","*.pdf","*.fai","*.mmi"} -e "ssh -p 22" --progress "../../../genomes_and_annotation/" "tomsrasp@147.156.206.251:/var/www/html/genomes_and_annotation/"

# echo transferring to VitViz

# rsync -avz --delete --exclude={"*.py","*.pptx","*.pdf","*.fai","*.mmi"} -e "ssh -p 4144" --progress "../../../genomes_and_annotation/" "tomsrasp2@147.156.206.186:/media/tomsrasp2/disk_a/genomes_and_annotation/"

# echo transferring to Antonio

# rsync -avz --delete --exclude={"*.py","*.pptx","*.pdf","*.fai","*.mmi"} -e "ssh -p 4144" --progress "../../../genomes_and_annotation/" "tomslab2@147.156.207.73:/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/"

# echo transferring to Iñigo

# rsync -avz --delete --exclude={"*.py","*.pptx","*.pdf","*.fai","*.mmi"} -e "ssh -p 4144" --progress "../../../genomes_and_annotation/" "inigo@147.156.206.222:/media/inigo/F4DE8D50DE8D0BD41/inigo/genomes_and_annotation/"

echo transferring to Garnatxita

rsync -avz --exclude={"*.py","*.pptx","*.pdf","*.fai","*.mmi"} -e "ssh -p 4144" --progress "../../../genomes_and_annotation/" "garnatxita@147.156.207.235:/media/hard_drive1/genomes_and_annotation/"
