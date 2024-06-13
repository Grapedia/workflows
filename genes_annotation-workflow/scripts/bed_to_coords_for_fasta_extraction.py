#!/usr/bin/env python
# coding: utf-8
import argparse


def command_line():
    parser = argparse.ArgumentParser(description='''
            Convert BED file to TXT file containing coordinates in the following format: chr:start-end
        '''.replace(12*' ', ''), formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='BED file', required=True)
    parser.add_argument('-o', '--output', help='Output file',
        required=True)
    return parser.parse_args()


def bed_to_coords(input_file, output):
    lines = []
    with open(input_file, 'r') as opened_file:
        for line in opened_file:
            lines.append('{0}:{1}-{2}'.format(*line[:-1].split('\t')))

    with open(output, 'w') as opened_file:
        opened_file.write('\n'.join(lines))


if __name__ == '__main__':
    args = command_line()
    bed_to_coords(args.input, args.output)
