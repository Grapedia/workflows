#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import re

if len(sys.argv) != 3:
    print("Usage: clean_protein_fasta.py input.fasta output.fasta")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

records = []
for i, record in enumerate(SeqIO.parse(input_file, "fasta"), start=1):
    # Clean the sequence ID
    simple_id = re.split(r"[ \t|]", record.id)[0]
    new_id = f"protein{i:06d}"
    
    # Clean AA sequence : replace '.' or '*' by 'X'
    cleaned_seq = str(record.seq).replace('.', 'X').replace('*', 'X').upper()

    # Create a new validated record
    clean_record = SeqRecord(
        seq=Seq(cleaned_seq),
        id=new_id,
        description=""
    )
    records.append(clean_record)

# Write the final cleaned file for BRAKER3
with open(output_file, "w") as out_f:
    SeqIO.write(records, out_f, "fasta")
