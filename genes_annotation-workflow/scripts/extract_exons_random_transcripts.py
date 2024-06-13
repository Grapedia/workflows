#!/usr/bin/env python
# coding: utf-8
import argparse
import gffutils
import random
import os


def command_line():
    parser = argparse.ArgumentParser(description='''
        Select n random transcripts from PsiCLASS transcriptome given in entry
        and extract their exons coodinates in TSV format.
        '''.replace(12*' ', ''), formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='PsiCLASS transcriptome in GFF3 format')
    parser.add_argument('-n', '--numberTranscripts', help='Number of transcript to extract')
    parser.add_argument('-o', '--output', help='Output file in TSV format')
    return parser.parse_args()


# Feature table creation based on the GFF3 given in input
def create_database(input_file):
    return gffutils.create_db(input_file, 'gff3_db', merge_strategy='create_unique', keep_order=True)


# Store the list of every transcript ID
def get_transcripts_id(db):
    return [f.id for f in db.features_of_type('transcript')]


# Get subset of transcripts ID based on the value given in input with the option "--numberTranscripts"
def get_random_transcript(list_id, n):
    indices = random.sample(range(len(list_id)), int(n))
    return [list_id[i] for i in sorted(indices)]


# Extract exons coordinates for each transcript given in input of this function
def extract_exons(db, list_id):
    list_by_id = []
    # For each transcript get their exons
    for i in list_id:
        list_children = []
        for c in db.children(i, order_by='start'):
            list_children.append('{0}\t{1}\t{2}'.format(c.seqid, c.start, c.end))
        list_by_id.append('\n'.join(list_children))
    return list_by_id


def write_output(list_lines, output_file):
    with open(output_file, 'w') as opened_file:
        opened_file.write('\n\n'.join(list_lines))


if __name__ == '__main__':
    args = command_line()
    gff3_db = create_database(args.input)
    list_transcripts_id = get_transcripts_id(gff3_db)
    if args.numberTranscripts:
        subset_transcripts_id = get_random_transcript(list_transcripts_id, args.numberTranscripts)
        list_coords = extract_exons(gff3_db, subset_transcripts_id)
    else:
        list_coords = extract_exons(gff3_db, list_transcripts_id)
    write_output(list_coords, args.output)
    os.remove('gff3_db')
