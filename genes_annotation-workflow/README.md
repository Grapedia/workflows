# TITAN: The Intensive Transcript Annotation Pipeline

📖 **Full documentation**: [Read the Docs](https://grapedia.readthedocs.io/en/latest/workflows.html#titan-the-intensive-transcript-annotation-pipeline)

## **Inputs**

The following parameters are defined in `nextflow.config` and are required for the pipeline execution.

### **containerOptions**
- **`containerOptions`**: Path to the **workdir pathway** (the directory where there is the main.nf script).  
  _Example_: `containerOptions = "-v /path/to/workflows/genes_annotation-workflow:/path/to/workflows/genes_annotation-workflow"`

### **Genome Assemblies**
- **`previous_assembly`**: Path to the **reference genome assembly** (FASTA).  
  _Example_: `data/assemblies/T2T_ref.fasta`
- **`new_assembly`**: Path to the **new genome assembly** to be annotated (FASTA).  
  _Example_: `data/assemblies/riesling.hap1.chromosomes.phased.fa`

### **Annotations & Data evidences**
- **`previous_annotations`**: Path to the **GFF3 file** containing annotations for the previous assembly.  
  _Example_: `data/annotations/PN40024_5.1_on_T2T_ref_with_names.gff3`
- **`RNAseq_samplesheet`**: Path to the **RNA-seq samplesheet** (CSV file) listing RNA-seq datasets to be used.  
  _Example_: `data/RNAseq_data/RNAseq_samplesheet.csv`
- **`protein_samplesheet`**: Path to the **protein data samplesheet** (CSV file) listing protein datasets to be used.  
  _Example_: `data/protein_data/samplesheet.csv`

### **Optional Parameters**
- **`EDTA`**: Whether to run **EDTA (transposable element annotation tool)**.  
  _Options_: `"yes"` or `"no"` (default: `"no"`)
  ⚠️ If EDTA is not set to "yes", the Aegis and Diamond2Go steps will not be executed, as Aegis requires a hard-masked genome. TITAN will only generate the evidence data.
- **`use_long_reads`**: Flag to indicate whether **long-read sequencing data** should be used.  
  _Options_: `true` or `false` (default: `true`)

## Workflow DAG

![Workflow Diagram](data_example/workflow.dag.jpg)

## Workflow Components

### **Liftoff Annotations**
Generates a **GFF3 file**, which is used by **Aegis** in the final step.

### **StringTie Merging (Short Reads - HISAT2)**
Generates **GTF file(s)**.  
- For each **short-read RNA-seq dataset**, a transcriptome is assembled using **HISAT2/StringTie**.  
- Transcriptomes are then **merged separately** for **stranded** and **unstranded** samples (unstranded is optional).  
- ⚠️ **This output is not used yet**—it is intended for **lncRNA detection**, which is not yet implemented.

### **StringTie Merging (Short Reads - STAR)**
Generates **GTF file(s)**, which are used by **Aegis** in the final step.  
- For each **short-read RNA-seq dataset**, a transcriptome is assembled using **STAR/StringTie** with **default** and **alternative (alt) parameters**.  
- Transcriptomes are **merged separately** for **stranded** and **unstranded** samples (unstranded is optional) and for **default and alt parameters** separately.

### **StringTie Merging (Long Reads)**
Generates a **single GTF file**, which is used by **Aegis** in the final step.  
- For each **IsoSeq RNA-seq dataset**, a transcriptome is assembled using **Minimap2/StringTie** with **default** and **alt parameters**.  
- These transcriptomes are then **merged into a single GTF file**, for **default and alt parameters** separately.
- ⚠️ **Long-read integration is optional**.

### **BRAKER3 Gene Prediction**
#### **braker3_prediction OR braker3_prediction_with_long_reads**
- **BRAKER3** is run with all **BAM files**, generated using **STAR (short reads)** or **Minimap2 (long reads)**.
- It also incorporates **protein FASTA files** to generate:
  - A **GeneMark GTF file**.
  - An **AUGUSTUS GFF3 file**, both of which are used by **Aegis** in the final step.

### **GFFCompare**
Generates a **GTF file**, which is used by **Aegis** in the final step.  
- For each **short-read RNA-seq dataset**, a transcriptome is assembled using **STAR/PsiCLASS**.  
- Transcriptomes are **merged separately** for **stranded** and **unstranded** samples using **GFFCompare** (unstranded is optional).

### **EDTA**
*(Integrated into TITAN, though not explicitly shown in the DAG)*  
Generates a **FASTA file** of the **hard-masked genome assembly**, which is used by **Aegis** in the final step.  

- **Input**: The **FASTA file** of the genome assembly to be annotated.  
- **Output**: A **hard-masked version** of the genome assembly.  

---

## **Final Integration by Aegis**
At the final step, **Aegis** integrates multiple annotation sources:  
✔ **Hard-masked genome assembly**
✔ **Liftoff annotations**  
✔ **GeneMark annotations**  
✔ **AUGUSTUS annotations**  
✔ **Stranded (and unstranded, if available) annotations from STAR/StringTie**  
✔ **IsoSeq annotations (if available)**  
✔ **Stranded (and unstranded, if available) annotations from STAR/PsiCLASS**  
✔ **protein FASTA files**

### **📖 Reference**  
Soon available ...

This results in the final **GFF3 annotation file**.

## **Diamond2GO**  
*(Integrated into TITAN, though not explicitly shown in the DAG)*  

Diamond2GO performs **functional gene annotation** based on the **final Aegis protein FASTA files** (both "main" and "all" sets).  

### **🔹 Input**  
- FASTA files containing the **"main"** and **"all"** protein sequences generated by Aegis.  

### **🔹 Output**  
- Functional annotations assigned by **Diamond2GO**.  

## **📖 Reference**  

**Aegis** : TO DO

**BRAKER3** : Gabriel et al. *BRAKER3: Fully automated genome annotation using RNA-seq and protein evidence with GeneMark-ETP, AUGUSTUS, and TSEBRA*  
📄 [Genome Res, 2024](https://pubmed.ncbi.nlm.nih.gov/38866550/)

**Diamond2GO** : Golden et al. *DIAMOND2GO: A rapid Gene Ontology assignment and enrichment tool for functional genomics.*  
📄 [bioRxiv, 2024](https://www.biorxiv.org/content/10.1101/2024.08.19.608700v1)

**EDTA** : Ou et al. *Benchmarking Transposable Element Annotation Methods for Creation of a Streamlined, Comprehensive Pipeline.*  
📄 [Genome Biol., 2019](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1905-y)

**fastp** : Chen et al. *fastp: an ultra-fast all-in-one FASTQ preprocessor*  
📄 [Bioinformatics., 2018](https://academic.oup.com/bioinformatics/article/34/17/i884/5093234)

**GFFCompare** : Pertea et al. *GFF Utilities: GffRead and GffCompare.*  
📄 [F1000Research., 2020](https://pubmed.ncbi.nlm.nih.gov/32489650/)

**GffRead** : Pertea et al. *GFF Utilities: GffRead and GffCompare.*  
📄 [F1000Research., 2020](https://pubmed.ncbi.nlm.nih.gov/32489650/)

**HISAT2** : Kim et al. *Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype*  
📄 [Nature biotechnology, 2019](https://www.nature.com/articles/s41587-019-0201-4)

**Liftoff** : Shumate et al. *Liftoff: accurate mapping of gene annotations*  
📄 [Bioinformatics, 2021](https://academic.oup.com/bioinformatics/article/37/12/1639/6035128?login=false)

**Minimap2** : Li et al. *Minimap2: pairwise alignment for nucleotide sequences*  
📄 [Bioinformatics, 2018](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778)

**PsiCLASS** : Song et al. *A multi-sample approach increases the accuracy of transcript assembly*  
📄 [Nature communications., 2019](https://www.nature.com/articles/s41467-019-12990-0)

**sra-tools** : [GitHub repository](https://github.com/ncbi/sra-tools)

**STAR** : Dobin et al. *STAR: ultrafast universal RNA-seq aligner*  
📄 [Bioinformatics, 2013](https://academic.oup.com/bioinformatics/article/29/1/15/272537)

**Salmon** : Patro et al. *Salmon provides fast and bias-aware quantification of transcript expression*  
📄 [Nature methods, 2017](https://www.nature.com/articles/nmeth.4197)

**Samtools** : Li et al. *The Sequence Alignment/Map format and SAMtools*  
📄 [Bioinformatics, 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2723002/)

**StringTie** : Shumate et al. *Improved transcriptome assembly using a hybrid of long and short reads with StringTie*  
📄 [PLOS Computational Biology, 2022](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009730)


