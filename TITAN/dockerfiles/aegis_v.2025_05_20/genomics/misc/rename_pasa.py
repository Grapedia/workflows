# Vitvi05_01chr07g00000_CS

input_gff = '/media/tomslab2/Storage/Antonio/TOMSBioLab/Camille_T2T/Rename_PASA/sample_mydb_pasa.sqlite.gene_structures_post_PASA_updates.1981355.gff3'
output_gff = '/media/tomslab2/Storage/Antonio/TOMSBioLab/Camille_T2T/Rename_PASA/test.gff3'

def extract_geneID(line):

    for field in line[-1].split(';'):

        if 'ID=' in field:

            geneID = field.replace('ID=', '')

            return geneID

def append_to_dictionary(dictionary, key, value):

    if key not in dictionary:

        dictionary[key] = [value]

    else:

        dictionary[key].append(value)

def parse_gff(input_gff):

    gff_genes = {}

    with open(input_gff, 'r') as f:

        for line in f:

            line = line.strip()
            if (line != '') and (line[0] != '#'):

                line_split = line.split('\t')

                if line_split[2] == 'gene':

                    geneID = extract_geneID(line_split)

                append_to_dictionary(gff_genes, geneID, line)

    return gff_genes


def change_geneIDs(gff_genes):

    new_gff_genes = {}
    gene_counter = 1

    for gene in gff_genes:

        mRNA_counter = 1
        exon_counter = 1
        CDS_counter = 1

        for line in gff_genes[gene]:

            line = line.split('\t')
            chromosome = line[0]
            structure = line[2]
            
            if structure == 'gene':

                geneID = f'Vitvi05_01{chromosome}g05{str(gene_counter).zfill(3)}'
                attribute = f'ID={geneID}'
                new_gff_genes[geneID] = []
                gene_counter += 1

            elif structure == 'mRNA':

                mRNAID = f'{geneID}.{mRNA_counter}'
                attribute = f'ID={mRNAID}; Parent={geneID}'
                mRNA_counter += 1

            elif structure == 'exon':

                exonID = f'{mRNAID}.exon{exon_counter}'
                attribute = f'ID={exonID}; Parent={mRNAID}'
                exon_counter += 1

            elif structure == 'CDS':

                CDSID = f'{mRNAID}.CDS{CDS_counter}'
                attribute = f'ID={CDSID}; Parent={mRNAID}'
                CDS_counter += 1

            elif structure == 'five_prime_UTR':

                continue

            elif structure == 'three_prime_UTR':

                continue

            else:

                print(f'Error: Unknown structure: "{structure}"')

            line[-1] = attribute
            line = '\t'.join(line)
            new_gff_genes[geneID].append(line)

    
    return new_gff_genes

def write_output_gff(new_gff_genes, output_path):

    with open(output_path, 'w') as o:

        for gene in new_gff_genes:

            for line in new_gff_genes[gene]:

                o.write(f'{line}\n')


gff_genes = parse_gff(input_gff)
new_gff_genes = change_geneIDs(gff_genes)
write_output_gff(new_gff_genes, output_gff)
                
