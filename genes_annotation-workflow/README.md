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
- **`use_long_reads`**: Flag to indicate whether **long-read sequencing data** should be used.  
  _Options_: `true` or `false` (default: `true`)

### **Logging**
- **`logfile`**: Path to the log file where pipeline execution logs will be saved.  
  _Example_: `pipeline_execution.log`

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

### **EDTA** (Not shown in the DAG but integrated into TITAN)
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

This results in the final **GFF3 annotation file**.
