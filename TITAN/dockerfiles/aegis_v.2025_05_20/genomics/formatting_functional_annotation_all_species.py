#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2024/05/28

Summarising functional annotation data.

@author: David Navarro
"""

import gc
from config.paths import all_species
from importlib import import_module
import os

subset_species = []

# redo summary files even if they already exist
overwrite = False

print(f"These are all the config files/species to deal with: {all_species}")

for config_species in all_species:
    if subset_species != []:
        if config_species not in subset_species:
            continue

    gc.collect()
    module = import_module(f"config.{config_species}")
    globals().update(vars(module))

    aegis_output = f"{path}aegis_output/"
    functional_annotation_path = f"{aegis_output}functional_annotation_output/"

    os.system(f"mkdir -p {functional_annotation_path}")

    print(f"Working on {config_species}")

    for annotation in annotation_files:
        annotation_found = False
        for file in functional_annotation:
            if annotation in file:
                annotation_found = True

        mapman_file = ""
        interpro_file = ""
        eggnog_file = ""
        gff_file = annotation_files[annotation]

        if annotation_found:
            for file in functional_annotation:
                if annotation in file:
                    if file.split("_")[-1] == "mapman":
                        mapman_file = functional_annotation[file]
                    elif file.split("_")[-1] == "eggnog":
                        eggnog_file = functional_annotation[file]
                    elif file.split("_")[-1] == "interpro":
                        interpro_file = functional_annotation[file]


            out_summary = f"{functional_annotation_path}{annotation}_functional_annotation_summary.csv"

            if not os.path.isfile(out_summary) or overwrite:

                genes = {}
                f_in = open(gff_file)
                for line in f_in:
                    if line.startswith("#"):
                        continue
                    line = line.strip().split("\t")
                    if len(line) < 3:
                        continue
                    elif line[2] == "gene" or line[2] == "pseudogene" or line[2] == "transposable_element_gene":
                        attributes = line[-1].split(";")
                        for a in attributes:
                            if "ID=" in a:
                                ID = a.split("=")[1].strip()
                                if ID not in genes:
                                    genes[ID] = { "eggNOG_descriptor" : set(), "mapman" : set(), "mapman_detailed" : set(), "Pfam" : set(), "PANTHER" : set(), "interpro" : set(), "other_domain_databases" : set(), "GO_molecular_function" : set(), "GO_biological_process" : set(), "GO_cellular_component" : set(), "EC" : set(), "COG_category" : set(), "KEGG_ko" : set(), "KEGG_pathway" : set(), "KEGG_module" : set(), "KEGG_reaction" : set(), "KEGG_rclass" : set(), "KEGG_BRITE" : set(), "KEGG_TC" : set()}
                f_in.close()

                mapman_terms = {}
                if mapman_file != "":
                    f_in = open(mapman_file)
                    for line in f_in:
                        line = line.strip().split("\t")
                        code = line[0]
                        name = line[1]
                        if code not in mapman_terms:
                            mapman_terms[code] = name
                        mapman_genes = line[2:]
                        for gene in mapman_genes:
                            if config_species == "mulberry" or config_species == "tomato":
                                gene = gene.split(".")[0]
                            elif config_species == "chlamydomonas":
                                if gene.startswith("Cremt.") or gene.startswith("Crecp."):
                                    gene = gene[:3] + gene[3].upper() + gene[4:]
                            elif config_species == "snapdragon":
                                if gene.startswith("am") or gene.startswith("Am"):
                                    gene = gene.split(".")[0]
                                    gene = gene[0].upper() + gene[1:]
                            if gene in genes:
                                if name.split(".")[-1] == "No Mercator4 annotation" or name.split(".")[-1] == "no other annotation available" or name.split(".")[-1] == "not annotated" or name.split(".")[-1] == "annotated" or name.split(".")[-1] == "not assigned":
                                    continue
                                genes[gene]["mapman"].add(code)
                            else:
                                print(f"Error: Mapman gene {gene} not found in gff")
                    f_in.close()

                for gene in genes:
                    genes[gene]["mapman"] = list(genes[gene]["mapman"])
                    genes[gene]["mapman_detailed"] = genes[gene]["mapman"].copy()

                # simplifying mapman codes:
                total_codes_deleted = 0
                for gene in genes:
                    all_codes = genes[gene]["mapman"]
                    codes_to_delete = set()
                    for c1 in all_codes:
                        for c2 in all_codes:
                            if c1 != c2:
                                if c1 in c2 and len(c1) < len(c2) and c1 == c2[:len(c1)]:
                                    codes_to_delete.add(c1)
                                elif c2 in c1 and len(c2) < len(c1) and c2 == c1[:len(c2)]:
                                    codes_to_delete.add(c2)
                    total_codes_deleted += len(codes_to_delete)
                    for c in codes_to_delete:
                        genes[gene]["mapman"].remove(c)

                print(f"A total of {total_codes_deleted} mapman codes were removed in the 'mapman_simplified' column.")

                for gene in genes:
                    new_terms = []
                    for mterm in genes[gene]["mapman"]:
                        name = mapman_terms[mterm]
                        shortname = name.split(".")[-1]
                        new_terms.append(f"{mterm}: '{shortname}'")
                    genes[gene]["mapman"] = new_terms.copy()

                    new_terms = []
                    for mterm in genes[gene]["mapman_detailed"]:
                        name = mapman_terms[mterm]
                        new_terms.append(f"{mterm}: '{name}'")
                    genes[gene]["mapman_detailed"] = new_terms.copy()

                # making an ID dictionary for interpro results
                categories = {}
                #unacceptable_categories = ["Coils", "MobiDBLite"]
                unacceptable_categories = []

                if interpro_file != "":
                    f_in = open(interpro_file)
                    for line in f_in:
                        line = line.strip().split("\t")
                        gene = line[0].split("_t")[0]
                        if config_species == "tomato":
                            gene = gene.split(".")[0]

                        elif config_species == "snapdragon":
                            if gene.startswith("am") or gene.startswith("Am"):
                                gene = gene.split(".")[0]
                                gene = gene[0].upper() + gene[1:]
                            elif gene.startswith("AnM") or gene.startswith("anm"):
                                gene = gene.split(".")[0] + "." + gene.split(".")[1]
                                gene = gene[0].upper() + gene[1] + gene[2].upper() + gene[3].upper() + gene[4:]

                        cat = line[3]
                        ID = line[4]
                        name1 = line[5]
                        interpro = line[11]
                        name2 = line[12]
                        score = line[8]

                        if cat not in unacceptable_categories:

                            genes[gene]["interpro"].add((ID, name1, name2, score, cat))

                            if cat not in categories:
                                categories[cat] = {ID:(name1, name2, interpro)}
                            else:
                                if ID not in categories[cat]:
                                    categories[cat][ID] = (name1, name2, interpro)
                                else:
                                    if (name1, name2, interpro) != categories[cat][ID]:
                                        print(f"Error: non unique ID: {ID}, from category {cat} in interpro results")
                                
                    f_in.close()

                for gene in genes:
                    genes[gene]["interpro"] = list(genes[gene]["interpro"])
                    unique = []
                    for ID, name1, name2, score, cat in genes[gene]["interpro"]:
                        if ID in unique:
                            print(f"Warning: same {ID} but with different score ({score}) given to {gene}")
                            continue
                            
                        else:
                            combined_name = ""
                            if name1 != "-" and name2 != "-":
                                combined_name = f"{name1}-{name2}"
                            elif name1 != "-":
                                combined_name = name1
                            elif name2 != "-":
                                combined_name = name2
                            if combined_name != "":
                                unique.append(f"{ID} ('{combined_name}', '{score}', '{cat}')")
                            else:
                                unique.append(f"{ID} ('-', '{score}', '{cat}')")

                    genes[gene]["interpro"] = unique

                    for u in genes[gene]["interpro"]:
                        category = u.split("', '")[-1][:-2]
                        if category == "Pfam":
                            genes[gene]["Pfam"].add(u)
                        elif category == "PANTHER":
                            genes[gene]["PANTHER"].add(u)
                        else:
                            genes[gene]["other_domain_databases"].add(u)

                    del genes[gene]["interpro"]

                go_terms = {}

                f_in = open(GO_file)
                pre_pre_previous_line = ""
                pre_previous_line = ""
                previous_line = ""
                for line in f_in:
                    line = line.strip()
                    if pre_pre_previous_line == "[Term]" and "id: " in pre_previous_line and "name: " in previous_line:
                        id = pre_previous_line.split(": ")[1]
                        name = previous_line.split(": ")[1]
                        namespace = line.split(": ")[1]
                        if id not in go_terms:
                            go_terms[id] = (name, namespace)
                        else:
                            print("Warning: repeat ids in go database file")
                    pre_pre_previous_line = pre_previous_line
                    pre_previous_line = previous_line
                    previous_line = line
                f_in.close()

                if eggnog_file != "":
                    f_in = open(eggnog_file)
                    for line in f_in:
                        if line.startswith("#"):
                            continue
                        line = line.strip().split("\t")

                        if config_species == "tomato":
                            id = line[0].split(".")[0]
                        
                        elif config_species == "snapdragon":
                            if line[0].startswith("am") or line[0].startswith("Am"):
                                id = line[0].split(".")[0]
                                id = id[0].upper() + id[1:]
                            elif line[0].startswith("AnM") or line[0].startswith("anm"):
                                id = line[0].split(".")[0] + "." + line[0].split(".")[1]
                                id = id[0].upper() + id[1] + id[2].upper() + id[3].upper() + id[4:]
                        else:
                            id = line[0].split("_t")[0]

                        egg_nog_descriptor = line[7]
                        egg_nog_symbol = line[8]

                        egg_nog = ""
                        if egg_nog_descriptor != "-":
                            egg_nog = egg_nog_descriptor
                        if egg_nog_symbol != "-":
                            egg_nog += f" ({egg_nog_symbol})"

                        GOs = set(line[9].split(","))

                        ECs = set(line[10].split(","))
                        ECs.discard("-")
                        COG_category = set(line[6].split(","))
                        COG_category.discard("-")
                        KEGG_ko = set(line[11].split(","))
                        KEGG_ko.discard("-")
                        KEGG_pathway = set(line[12].split(","))
                        KEGG_pathway.discard("-")
                        KEGG_module = set(line[13].split(","))
                        KEGG_module.discard("-")
                        KEGG_reaction = set(line[14].split(","))
                        KEGG_reaction.discard("-")
                        KEGG_rclass = set(line[15].split(","))
                        KEGG_rclass.discard("-")
                        KEGG_BRITE = set(line[16].split(","))
                        KEGG_BRITE.discard("-")
                        KEGG_TC = set(line[17].split(","))
                        KEGG_TC.discard("-")

                        for GO in GOs:
                            if GO in go_terms:
                                go_class = go_terms[GO][1]
                                if go_class == "biological_process":
                                    genes[id]["GO_biological_process"].add(f"{GO}: '{go_terms[GO][0]} ({go_class})'")
                                elif go_class == "cellular_component":
                                    genes[id]["GO_cellular_component"].add(f"{GO}: '{go_terms[GO][0]} ({go_class})'")
                                elif go_class == "molecular_function":
                                    genes[id]["GO_molecular_function"].add(f"{GO}: '{go_terms[GO][0]} ({go_class})'")
                                else:
                                    print(go_class)

                        if id in genes:
                            genes[id]["eggNOG_descriptor"].add(egg_nog)
                            genes[id]["EC"] = ECs
                            genes[id]["COG_category"] = COG_category
                            genes[id]["KEGG_ko"] = KEGG_ko
                            genes[id]["KEGG_pathway"] = KEGG_pathway
                            genes[id]["KEGG_module"] = KEGG_module
                            genes[id]["KEGG_reaction"] = KEGG_reaction
                            genes[id]["KEGG_rclass"] = KEGG_rclass
                            genes[id]["KEGG_BRITE"] = KEGG_BRITE
                            genes[id]["KEGG_TC"] = KEGG_TC

                    f_in.close()

                out = ""

                for gene in genes:
                    header = f"gene"
                    for key in genes[gene]:
                        header += f"\t{key}"
                    break
                out = f"{header}\n"

                for gene in genes:
                    outline = gene
                    for key in genes[gene]:
                        values = "; ".join(genes[gene][key])
                        outline += f"\t{values}"
                    out += f"{outline}\n"

                f_out = open(out_summary, "w", encoding="utf-8")
                f_out.write(out)
                f_out.close()

