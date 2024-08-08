# Genes annotation workflow

```
----- INPUT
previous_assembly : <FASTA_FILE>
new_assembly : <FASTA_FILE>
previous_annotations : <GFF3_FILE>
RNAseq_samplesheet : <CSV_FILE>
protein_samplesheet : <CSV_FILE>
---- OUTPUT
output file : <GFF3_FILE>
```
## Input data structure :

```
├── data
│   ├── annotations
│   ├── assemblies
│   ├── protein_data
│   └── RNAseq_data
```

data/annotations : contains the previous annotation in GFF3 format (eg : Vitis_vinifera_gene_annotation_on_V2_20.gff3)

data/assemblies : contains previous assembly (eg PN12Xv2.fasta) and new assembly (eg Chinese_ref_v2.fa)

data/protein_data : contains all the protein data files (FASTA) to perform protein alignments. Contains also a samplesheet describing the protein data file to use.

          Example :

          ```
          organism,filename,maker_braker2
          arabidopsis,arabidopsis_prot_2022_01.fasta,no
          viridiplantae,Viridiplantae_swissprot.fasta,yes
          eudicotyledones_uniprot,eudicotyledons_uniprot.fasta,no
          eudicotyledones_orthoDB,eudicotyledons_odb10.fasta,yes
          vitales,vitales.fasta,no
          ```

data/RNAseq_data : contains all the RNAseq data for transcriptome assembly. Contains also the RNAseq_samplesheet. If FASTQ, the fastq file must be in the right folder, if SRA, the workflow will download the SRA file and convert it to fastq.gz file.

          Example of RNAseq_samplesheet :

          ```
          sampleID,SRA_or_FASTQ,paired_or_single
          ERR1059552,FASTQ,paired
          ERR1059553,FASTQ,paired
          ERR1059554,SRA,paired
          ERR1059555,SRA,paired
          SRR5435969,FASTQ,paired
          SRR8775072,FASTQ,paired
          SRR3046429,SRA,paired
          SRR3046438,SRA,paired
          SRR520373,SRA,single
          ```

WARNING : in data/RNAseq_data, for the FASTQ files, the name need to be ${sampleID}.fastq.gz for single-end and ${sampleID}_1.fastq.gz and ${sampleID}_2.fastq.gz for paired-end.
