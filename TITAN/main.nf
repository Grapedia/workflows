nextflow.enable.dsl = 2

params.RNAseq_samplesheet = false
params.protein_samplesheet = false
params.new_assembly = false
params.previous_assembly = false
params.previous_annotations = false
params.output_dir = false

if (!params.output_dir || !params.RNAseq_samplesheet || !params.protein_samplesheet || !params.new_assembly || !params.previous_assembly || !params.previous_annotations) {
    error "Missing required parameters. Please provide values for: output_dir, RNAseq_samplesheet, protein_samplesheet, new_assembly, previous_assembly, previous_annotations"
}

// Include subworkflows
include { generate_evidence_data } from './subworkflows/generate_evidence_data'
include { aegis } from './subworkflows/aegis'

// Define the output directory for intermediate files
params.workflow = params.workflow ?: "all"  // default = "all"

// Define Channels
Channel.fromPath(file(params.RNAseq_samplesheet))
       .splitCsv(header: true, sep: ',')
       .filter( ~/.*long.*/ )
       .map { row -> [ row.sampleID, row.SRA_or_FASTQ, row.library_layout ] }
       .set{ samples_list_long_reads }


Channel.fromPath(file(params.RNAseq_samplesheet))
       .splitCsv(header: true, sep: ',')
       .filter( ~/.*single.*/ )
       .map { row -> [ row.sampleID, row.SRA_or_FASTQ, row.library_layout ] }
       .set{ samples_list_single_short_reads }


Channel.fromPath(file(params.RNAseq_samplesheet))
       .splitCsv(header: true, sep: ',')
       .filter( ~/.*paired.*/ )
       .map { row -> [ row.sampleID, row.SRA_or_FASTQ, row.library_layout ] }
       .set{ samples_list_paired_short_reads }

// Combined single and paired-end reads
samples_list_single_short_reads
    .concat(samples_list_paired_short_reads)
    .set { samples_list_short_reads }

Channel.fromPath(file(params.protein_samplesheet))
       .splitCsv(header: true, sep: ',')
       .map { row -> [ row.organism, row.filename ] }
       .set{ protein_list }

workflow {
    // if the workflow parameter is equal to generate_evidence_data or all, the first workflow is executed
    if (params.workflow == "generate_evidence_data" || params.workflow == "all") {
        evidence_data = generate_evidence_data(
            samples_list_long_reads, 
            samples_list_short_reads, 
            protein_list
        )
        
        if (params.workflow == "all") {
            aegis(evidence_data.results)
        }
    }
    
    // if the workflow parameter is equal to aegis, the outputs generated by the first workflow (generate_evidence_data) must be retrieved manually.
    else if (params.workflow == "aegis") {

      def outdir_1 = params.output_dir
      println "Checking output directory: ${outdir_1}"

      def input_files = file("${outdir_1}").listFiles()

      // if (!input_files || input_files.size() == 0) {
      //     println "ERROR : No file found in ${outdir_1} !"
      //     exit(1)
      // } else {
      //     input_files.each { file ->
      //         println "Files found : ${file}"
      //     }
      // }

      def workflow_inputs = []

      // Load manually the outputs of generate_evidence_data workflow
      if (params.EDTA == 'yes') {
        if (file("${outdir_1}/assembly_masked.EDTA.fasta").exists()) {
          workflow_inputs << tuple("masked_genome.masked_genome", file("${outdir_1}/assembly_masked.EDTA.fasta"))
        } else {
          println "ERROR : File ${outdir_1}/assembly_masked.EDTA.fasta doesn't exists !"
        }
      } else {
        println "EDTA = No, we need this output for Aegis. EXIT."
      }

      if (params.use_long_reads) {
        if (file("${outdir_1}/merged_minimap2_stringtie_long_reads_default.gtf").exists()) {
          workflow_inputs << tuple("merged_long_reads.default_args_gff", file("${outdir_1}/merged_minimap2_stringtie_long_reads_default.gtf"))
        } else {
          println "ERROR : ${outdir_1}/merged_minimap2_stringtie_long_reads_default.gtf doesn't exists !"
        }

        if (file("${outdir_1}/merged_minimap2_stringtie_long_reads_alt.gtf").exists()) {
          workflow_inputs << tuple("merged_long_reads.alt_args_gff", file("${outdir_1}/merged_minimap2_stringtie_long_reads_alt.gtf"))
        } else {
          println "ERROR : ${outdir_1}/merged_minimap2_stringtie_long_reads_alt.gtf doesn't exists !"
        }
      }

      if (file("${outdir_1}/augustus.hints.gff3").exists()) {
          workflow_inputs << tuple("braker3_results.augustus_gff", file("${outdir_1}/augustus.hints.gff3"))
      } else {
        println "ERROR : ${outdir_1}/augustus.hints.gff3 doesn't exists !"
      }
      if (file("${outdir_1}/genemark.gtf").exists()) {
          workflow_inputs << tuple("braker3_results.genemark_gtf", file("${outdir_1}/genemark.gtf"))
      } else {
        println "ERROR : ${outdir_1}/genemark.gtf doesn't exists !"
      }
      if (file("${outdir_1}/liftoff_previous_annotations.gff3").exists()) {
          workflow_inputs << tuple("previous_annotations.liftoff_previous_annotations", file("${outdir_1}/liftoff_previous_annotations.gff3"))
      } else {
        println "ERROR : ${outdir_1}/liftoff_previous_annotations.gff3 doesn't exists !"
      }
      if (file("${outdir_1}/merged_star_stringtie_stranded_default.gtf").exists()) {
          workflow_inputs << tuple("merged_star_stringtie.default_args_stranded", file("${outdir_1}/merged_star_stringtie_stranded_default.gtf"))
      } else {
        println "ERROR : ${outdir_1}/merged_star_stringtie_stranded_default.gtf doesn't exists !"
      }
      if (file("${outdir_1}/merged_star_stringtie_stranded_alt.gtf").exists()) {
          workflow_inputs << tuple("merged_star_stringtie.alt_args_stranded", file("${outdir_1}/merged_star_stringtie_stranded_alt.gtf"))
      } else {
        println "ERROR : ${outdir_1}/merged_star_stringtie_stranded_alt.gtf doesn't exists !"
      }
      if (file("${outdir_1}/merged_star_psiclass_stranded.gtf").exists()) {
          workflow_inputs << tuple("gffcompare_out.star_psiclass_stranded", file("${outdir_1}/merged_star_psiclass_stranded.gtf"))
      } else {
        println "ERROR : ${outdir_1}/merged_star_psiclass_stranded.gtf doesn't exists !"
      }
      if (file("${outdir_1}/merged_star_stringtie_unstranded_default.gtf").exists()) {
          workflow_inputs << tuple("merged_star_stringtie.default_args_unstranded", file("${outdir_1}/merged_star_stringtie_unstranded_default.gtf"))
      } else {
        println "WARNING : no RNAseq unstranded data used in this run !"
        tuple("merged_star_stringtie.default_args_unstranded", null)
      }
      if (file("${outdir_1}/merged_star_stringtie_unstranded_alt.gtf").exists()) {
          workflow_inputs << tuple("merged_star_stringtie.alt_args_unstranded", file("${outdir_1}/merged_star_stringtie_unstranded_alt.gtf"))
      } else {
        println "WARNING : no RNAseq unstranded data used in this run !"
        tuple("merged_star_stringtie.alt_args_unstranded", null)
      }
      if (file("${outdir_1}/merged_star_psiclass_unstranded.gtf").exists()) {
          workflow_inputs << tuple("gffcompare_out.star_psiclass_unstranded", file("${outdir_1}/merged_star_psiclass_unstranded.gtf"))
      } else {
        println "WARNING : no RNAseq unstranded data used in this run !"
        tuple("gffcompare_out.star_psiclass_unstranded", null)
      }

      def workflow_inputs_list = workflow_inputs // Keep the list as it is

      println "Files sent to Aegis : ${workflow_inputs_list}"

      aegis(workflow_inputs_list) // We send the list directly, not a Channel
  }
}
