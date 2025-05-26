#!/bin/bash


species=("V3Vitis" "Vesca")


for ((i = 0; i < ${#species[@]}; i++)); do
    for ((j = i + 1; j < ${#species[@]}; j++)); do
        echo "${species[i]} ${species[j]}"
        
        sed 's/_t[0-9][0-9][0-9]_CDS[0-9]*//g' raw_data/${species[i]}.fasta > output/${species[i]}_raw.cds
        sed 's/_t[0-9][0-9][0-9]_CDS[0-9]*//g' raw_data/${species[j]}.fasta > output/${species[j]}_raw.cds
        
        docker run -it --rm -v /media:/media -w $PWD jcvi-image2 python -m jcvi.formats.fasta format output/${species[i]}_raw.cds output/${species[i]}.cds
        docker run -it --rm -v /media:/media -w $PWD jcvi-image2 python -m jcvi.formats.fasta format output/${species[j]}_raw.cds output/${species[j]}.cds
        
        rm output/${species[i]}_raw.cds
        rm output/${species[j]}_raw.cds
        
        docker run -it --rm -v /media:/media -w $PWD jcvi-image2 python -m jcvi.formats.gff bed --type=mRNA --key=featurecounts_id --primary_only raw_data/${species[i]}.gff3 -o output/${species[i]}.bed
        docker run -it --rm -v /media:/media -w $PWD jcvi-image2 python -m jcvi.formats.gff bed --type=mRNA --key=featurecounts_id --primary_only raw_data/${species[j]}.gff3 -o output/${species[j]}.bed
        
        cd output
        docker run -it --rm -v /media:/media -w $PWD jcvi-image2 python -m jcvi.compara.catalog ortholog ${species[i]} ${species[j]} --no_strip_names
        cd ..
        
    done
done



# sed -i 's/_t[0-9][0-9][0-9]_CDS[0-9]*//g' ITAG4_all_CDSs.fasta
# python -m jcvi.formats.gff bed --type=mRNA --key=featurecounts_id --primary_only Vvinifera_145_Genoscope.12X.gene.gff3 -o grape.bed
# python -m jcvi.formats.fasta format Ppersica_298_v2.1.cds.fa.gz peach.cds
# docker run -it --rm -v /media:/media -w $PWD jcvi-image2 python -m jcvi.compara.catalog ortholog v51 ITAG4 --no_strip_names


