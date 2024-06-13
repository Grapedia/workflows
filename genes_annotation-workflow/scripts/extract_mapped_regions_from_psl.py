#!/usr/bin/env python
# coding: utf-8
import argparse


def command_line():
    parser = argparse.ArgumentParser(description='''
        Extract mapped regions on target assembly from alignments in blast8
        (tabular) format
        '''.replace(12*' ', ''), formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='Input file in BLAST8 format')
    parser.add_argument('-o', '--output', help='Output file in BED format')
    return parser.parse_args()


# Convert line in BED format: seqid  start  end
def convert_lines(line):
    full_list = line.split('\t')
    new_list = [full_list[i] for i in [13, 15, 16]]

    # If start > end -> invert
    if int(new_list[1]) > int(new_list[2]):
        val = new_list[1]
        new_list[1] = new_list[2]
        new_list[2] = val

    return '{0}\t{1}\t{2}'.format(*[l for l in new_list])


def parse_psl_file(input_file):
    mapped_regions = []
    with open(input_file, 'r') as opened_file:
        for align in opened_file:
            mapped_regions.append(convert_lines(align))
    return mapped_regions


def write_output(output_file, output_list):
    with open(output_file, 'w') as opened_file:
        opened_file.write('\n'.join(output_list))


if __name__ == '__main__':
    args = command_line()
    coords_mapped_regions = parse_psl_file(args.input)
    write_output(args.output, coords_mapped_regions)
