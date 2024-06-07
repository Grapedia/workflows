#!/usr/bin/env python3

# import packages
import argparse
import warnings
import gffutils
import os
from itertools import chain

def command_line():
  parser = argparse.ArgumentParser(description='''
      Create/write ID to each exon in gff3 output from PSIClass.
      ''')
  parser.add_argument('-g', '--gff3', help='GFF3 file from PSIClass.')
  parser.add_argument('-o', '--output', help='GFF3 final output file')
  args = parser.parse_args()
  return args.gff3, args.output


# creation of a gffutils database
def create_database(fgff3, db_filename):
  print("------------------")
  print("Creating", db_filename ," database...")
  gffutils.create_db(fgff3, db_filename, sort_attribute_values=True,merge_strategy="create_unique", force=True)
  print("------------------")


def create_and_add_exons_IDs(db_filename, finalgff3):
  output_file = open(finalgff3, 'w')
  db = gffutils.FeatureDB(db_filename)
  transcripts_iterator=chain(db.features_of_type('transcript'))
  for transcript in db.features_of_type('transcript'):
    output_file.write("{}\n".format(str(transcript)))
    count_exon=1
    #print(transcript.attributes['ID'][0])
    for exon in db.children(transcript.id, level=1):
      exon.attributes['ID']=str(transcript.attributes['ID'][0]) + ".exon" + str(count_exon)
      count_exon+=1
      output_file.write("{}\n".format(str(exon)))
  output_file.close()

if __name__ == '__main__':
  fgff3, foutput = command_line()

  # 1st step : create the databases from the Apollo's gff3 and PN40024.v4 ref
  create_database(fgff3, "db_genes_psiclass")
  create_and_add_exons_IDs("db_genes_psiclass",foutput)
