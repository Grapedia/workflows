from modules.geneclasses import Genome, Annotation

g = Genome("genome_tag", "genome_file_path")

a = Annotation("annotation_tag", "annotation_file_path", g)

# generating protein sequences
a.generate_sequences(g, just_CDSs=True)

# exporting main proteins
a.export_proteins(only_main=True, verbose=False, custom_path="output_folder")

# detecting gene overlaps
a_other = Annotation("2nd_annotation_tag", "2nd_annotation_file_path", g)
a.detect_gene_overlaps(a_other)
# overlaps table with scores
a.export_equivalences(stringent=False, return_df=False, NAs=False, export_csv=True, custom_path="output_folder")

masked_g = Genome("masked_genome_tag", "mask_genome_path")

# calculate masked fraction for each gene:

a.calculate_transcript_masking(hard_masked_genome=masked_g)

# at the moment we don't have an export method but this loop should print the masked_fraction property
# what is considered a masked gene is up to debate, at the moment we were only considering
# completely masked genes i.e. masked_fraction = 1

for genes in a.chrs.values():
    for g in genes.values():
        print(g.masked_fraction, g.id)