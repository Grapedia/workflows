#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from config.paths import root, GO_file

species="cannabis"
# assembly/annotation folder if any
assembly=""
species_path=f"{root}{species}/"
path = f"{species_path}{assembly}"

genome_files = {}
genome_files["CBDRx"] = f"{path}CBDRx_genome/CBDRx-18.genome.fasta"
genome_files["PurpleKush"] = f"{path}PurpleKush_genome/PurpleKush.genome.fasta"
genome_files["JamaicanLion"] = f"{path}jamaican_lion_genome/JamaicanLionDASH.genome.fasta"
genome_files["Finola"] = f"{path}finola_genome/Finola.genome.fasta"

# notice _on_ nomenclature
annotation_files = {}
annotation_files["CBDRx_on_CBDRx"] = f"{path}CBDRx_annotation/CBDRx-18_tidy.gff3"
annotation_files["PurpleKush_on_PurpleKush"] = f"{path}PurpleKush_annotation/PurpleKush_tidy.gff3"
annotation_files["JamaicanLion_on_JamaicanLion"] = f"{path}jamaican_lion_annotation/JamaicanLionDASH.gene_tidy.gff3"
annotation_files["Finola_on_Finola"] = f"{path}finola_annotation/Finola.gene.gff3"

# notice _from_ and _on_ nomenclature
# leave empty if no liftoff or equivalent files
annotation_transfer_files = {}

functional_annotation = {}
functional_annotation["JamaicanLion_on_JamaicanLion_mapman"] = f"{path}/jamaican_lion_annotation/functional_annotation/MapMan.Mercator.Jamaican.gmt"
functional_annotation["JamaicanLion_on_JamaicanLion_interpro"] = f"{path}/jamaican_lion_annotation/functional_annotation/JamaicanLion_aegis_on_JamaicanLion_all_proteins.interpro"
functional_annotation["JamaicanLion_on_JamaicanLion_eggnog"] = f"{path}/jamaican_lion_annotation/functional_annotation/out.emapper.annotations"

functional_annotation["PurpleKush_on_PurpleKush_mapman"] = f"{path}/PurpleKush_annotation/functional_annotation/MapMan.Mercator.PurpleKush.gmt"
functional_annotation["PurpleKush_on_PurpleKush_interpro"] = f"{path}/PurpleKush_annotation/functional_annotation/PurpleKush_aegis_on_PurpleKush_all_proteins.interpro"
functional_annotation["PurpleKush_on_PurpleKush_eggnog"] = f"{path}/PurpleKush_annotation/functional_annotation/out.emapper.annotations"

functional_annotation["CBDRx_on_CBDRx_mapman"] = f"{path}/CBDRx_annotation/functional_annotation/MapMan.Mercator.CBDRX.gmt"
functional_annotation["CBDRx_on_CBDRx_interpro"] = f"{path}/CBDRx_annotation/functional_annotation/CBDRx_aegis_on_CBDRx_all_proteins.interpro"
functional_annotation["CBDRx_on_CBDRx_eggnog"] = f"{path}/CBDRx_annotation/functional_annotation/out.emapper.annotations"

genome_dictionary = {}

for tag in genome_files:
    genome_dictionary[tag] = {}