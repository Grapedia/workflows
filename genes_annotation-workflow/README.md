# Genes annotation workflow

```
----- INPUT
assemblies_folder : <PATH>
previous_assembly : <FASTA_FILE>
new_assembly : <FASTA_FILE>
previous_annotations : <GFF3_FILE>
RNAseq_samplesheet : <CSV_FILE>
protein_samplesheet : <CSV_FILE>
geneid_param_file : <TAB_DELIMITED_FILE>
pasa_config_file : <TAB_DELIMITED_FILE>
---- OUTPUT
output file : <GFF3_FILE>
```
## Input data structure :

```
├── data
│   ├── annotations
│   ├── assemblies
│   ├── evidencemodeler_weights_file
│   ├── geneid_param_file
│   ├── pasa_config_file
│   ├── protein_data
│   └── RNAseq_data
│       ├── stranded
│       └── unstranded
```

data/annotations : contains the previous annotation in GFF3 format (eg : Vitis_vinifera_gene_annotation_on_V2_20.gff3)

data/assemblies : contains previous assembly (eg PN12Xv2.fasta) and new assembly (eg Chinese_ref_v2.fa)

data/evidencemodeler_weights_file : contains the weights.txt file for EvidenceModeler (to remove ?)

          Example :

          ```
          ABINITIO_PREDICTION maker 3
          ABINITIO_PREDICTION geneid_v1.4 1
          ABINITIO_PREDICTION GlimmerHMM 1
          ABINITIO_PREDICTION AUGUSTUS 1
          ABINITIO_PREDICTION Liftoff 5
          PROTEIN exonerate 1
          TRANSCRIPT  PsiCLASS_RNAseq_stranded  10
          TRANSCRIPT  PsiCLASS_RNAseq_unstranded  8
          ```

data/geneid_param_file : contains the parameter file for geneid (vvinifera.param.Jan_12_2007)

data/pasa_config_file : contains the config file for PASA (pasa.alignAssembly.Template.txt)

          Example :

          ```
          ## templated variables to be replaced exist as <__var_name__>

          # database settings
          DATABASE=pasa

          #######################################################
          # Parameters to specify to specific scripts in pipeline
          # create a key = "script_name" + ":" + "parameter"
          # assign a value as done above.

          #script validate_alignments_in_db.dbi
          validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=<__MIN_PERCENT_ALIGNED__>
          validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=<__MIN_AVG_PER_ID__>

          #script subcluster_builder.dbi
          subcluster_builder.dbi:-m=50
          ```

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
data/RNAseq_data/{stranded,unstranded} : contains all the RNAseq data for transcriptome assembly. Contains also the RNAseq_samplesheet. If FASTQ, the fastq file must be in the right folder, if SRA, the workflow will download the SRA file and convert it to fastq.gz file.

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

WARNING : in data/RNAseq_data/{stranded,unstranded}, for the FASTQ files, the name need to be ${sampleID}.fastq.gz for single-end and ${sampleID}_1.fastq.gz and ${sampleID}_2.fastq.gz for paired-end.
