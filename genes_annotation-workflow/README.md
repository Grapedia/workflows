# Genes annotation workflow

```
----- INPUT
assemblies_folder : <PATH>
previous_assembly : <FASTA_FILE>
new_assembly : <FASTA_FILE>
previous_annotations : <GFF3_FILE>
RNAseq_samplesheet : <CSV_FILE>
arabidopsis : <FASTA_FILE>
viridiplantae : <FASTA_FILE>
eudicotyledones_uniprot : <FASTA_FILE>
eudicotyledones_orthoDB : <FASTA_FILE>
vitales : <FASTA_FILE>
geneid_param_file : <TAB_DELIMITED_FILE>
pasa_config_file : <TAB_DELIMITED_FILE>
---- OUTPUT
output file : <GFF3_FILE>
```
Example of RNAseq_samplesheet :

```
sampleID,stranded_or_unstranded,SRA_or_FASTQ,paired_or_single
ERR1059552,stranded,FASTQ,paired
ERR1059553,stranded,FASTQ,paired
ERR1059554,stranded,SRA,paired
ERR1059555,stranded,SRA,paired
SRR5435969,unstranded,FASTQ,paired
SRR8775072,unstranded,FASTQ,paired
SRR3046429,unstranded,SRA,paired
SRR3046438,unstranded,SRA,paired
SRR520373,unstranded,SRA,single
```
