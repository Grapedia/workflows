#!/usr/bin/env python
# coding: utf-8
import argparse
import re


def command_line():
    parser = argparse.ArgumentParser(description='''
            Filter alignments by:
                * identity percentage (AveragePercentIdentity) > 25%
                * similarity percentage (AveragePercentSimilarity) > 50%
                * sequence alignment coverage (Alignment*100/QueryLen) > 80%
        '''.replace(12*' ', ''), formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', required=True, help='Exonerate GFF-like file')
    parser.add_argument('-o', '--output', required=True, help='GFF file filtered')
    return parser.parse_args()


def format_and_recalculate(gff_line):
    gff_line_list = gff_line.split('\t')
    if len(gff_line_list) == 9:
        # Extract start from mapped region id
        start_mapped_region = int(re.split(':|-', gff_line_list[0])[1])

        # Recalculate start and end of alignement to have the coordinates on the entire ...
        # ... chromosome and not only on the mapped region
        gff_line_list[3] = str(int(gff_line_list[3])+start_mapped_region-1)
        gff_line_list[4] = str(int(gff_line_list[4])+start_mapped_region-1)
    return '\t'.join(gff_line_list)


# Filter alignments based on identity, similarity and coverage
def filter_alignments(list_alignments, lines, scores):
    if scores['AveragePercentIdentity'] > 25 and \
    scores['AveragePercentSimilarity'] > 50 and \
    scores['Alignment']*100/scores['QueryLen'] > 80:
        [list_alignments.append(format_and_recalculate(l)) for l in lines]


def filter_and_convert(input_file, output_file):
    with open(input_file, 'r') as open_input_file:
        list_alignments = []
        list_gff_lines = []
        scores_dict = {
            'AveragePercentIdentity': 0,
            'AveragePercentSimilarity': 0,
            'Score': 0,
            'QueryLen': 0,
            'Alignment': 0
        }
        # For each block (which correspond to 1 alignment), we filter and extract ...
        # ... only the GFF lines
        is_gff_line = False
        list_gff_lines = []
        for line in open_input_file.readlines():
            if '# --- START OF GFF DUMP ---' in line:
                is_gff_line = True
                if len(list_gff_lines) > 0:
                    # Filter and append or not the previous alignment
                    filter_alignments(list_alignments, list_gff_lines, scores_dict)
                # Start list for the new alignment
                list_gff_lines = []
                list_gff_lines.append(line[:-1])

            # Means it's the end of GFF line, after this line are written the different ...
            # ... score values, and alignment and query length
            elif '# --- END OF GFF DUMP ---' in line:
                is_gff_line = False
                list_gff_lines.append(line[:-1])
            else:
                if not line.startswith(('Command line', 'Hostname')):
                    # We need to keep the lines beginning with an '#', otherwise ...
                    # ... the conversion in gff3 format won't work properly
                    if not line.startswith('#') and not is_gff_line:
                        score = line.split(':')
                        if len(score) > 1:
                            scores_dict[score[0]] = float(score[1][:-1])
                    list_gff_lines.append(line[:-1])

        # To filter the last alignment we have to call one last time the function ...
        # ... after every lines are parsed
        filter_alignments(list_alignments, list_gff_lines, scores_dict)

        # Write output file
        with open(output_file, 'w') as open_output_file:
            open_output_file.write('{}\n'.format('\n'.join(list_alignments)))


if __name__ == '__main__':
    args = command_line()
    filter_and_convert(args.input, args.output)
