nextflow.enable.dsl = 2

params.RNAseq_samplesheet = false
params.protein_samplesheet = false
params.new_assembly = false
params.previous_assembly = false
params.previous_annotations = false
params.output_dir = false
params.egapx_paramfile = false

// Include subworkflows
include { generate_evidence_data } from './subworkflows/generate_evidence_data'
include { aegis } from './subworkflows/aegis'

// Define the output directory for intermediate files
params.workflow = params.workflow ?: "generate_evidence_data"  // default = "generate_evidence_data"

def isMissingParam(value) {
    return value == null || value == true || value == false || value.toString().trim() == '' || value.toString().trim().equalsIgnoreCase('true') || value.toString().trim().equalsIgnoreCase('false')
}

def validateRequiredParams(requiredParams) {
    def missingParams = requiredParams.findAll { paramName -> isMissingParam(params[paramName]) }
    if (missingParams) {
        def formattedParams = missingParams.collect { "--${it}" }.join(', ')
        error "Missing required parameter(s): ${formattedParams}"
    }
}

def validateExistingInputFiles(requiredFiles) {
    def missingFiles = requiredFiles.findAll { paramName -> !file(params[paramName]).exists() }
    if (missingFiles) {
        def formattedFiles = missingFiles.collect { paramName -> "--${paramName}: ${params[paramName]}" }.join('\n  ')
        error "Required input file(s) not found:\n  ${formattedFiles}"
    }
}

def validateWorkflowName() {
    def allowedWorkflows = ['generate_evidence_data', 'aegis']
    if (!allowedWorkflows.contains(params.workflow)) {
        error "Invalid --workflow '${params.workflow}'. Allowed values: ${allowedWorkflows.join(', ')}"
    }
}

def normalizeBooleanParam(value, paramName) {
    if (value == null) {
        return false
    }
    if (value instanceof Boolean) {
        return value
    }

    def normalizedValue = value.toString().trim().toLowerCase()
    if (['true', 'yes', 'y', '1'].contains(normalizedValue)) {
        return true
    }
    if (['false', 'no', 'n', '0'].contains(normalizedValue)) {
        return false
    }

    error "Invalid boolean value for --${paramName}: '${value}'. Allowed values: true/false, yes/no, 1/0"
}

def normalizeRuntimeFlags() {
    def useLongReads = normalizeBooleanParam(params.use_long_reads, 'use_long_reads')
    def runEgapx = normalizeBooleanParam(params.run_egapx, 'run_egapx')
    def runEdta

    if (params.run_edta != null) {
        runEdta = normalizeBooleanParam(params.run_edta, 'run_edta')
    } else {
        runEdta = normalizeBooleanParam(params.EDTA, 'EDTA')
    }

    // Keep the historical EDTA parameter synchronized while downstream modules migrate.
    params.EDTA = runEdta ? 'yes' : 'no'

    return [
        run_edta: runEdta,
        run_egapx: runEgapx,
        use_long_reads: useLongReads
    ]
}

workflow {
    validateWorkflowName()
    validateRequiredParams(['output_dir', 'egapx_paramfile', 'RNAseq_samplesheet', 'protein_samplesheet', 'new_assembly', 'previous_assembly', 'previous_annotations'])
    validateExistingInputFiles(['egapx_paramfile', 'RNAseq_samplesheet', 'protein_samplesheet', 'new_assembly', 'previous_assembly', 'previous_annotations'])
    def runtime_flags = normalizeRuntimeFlags()

    // if the workflow parameter is equal to generate_evidence_data the first workflow is executed
    if (params.workflow == "generate_evidence_data") {
        Channel.fromPath(params.RNAseq_samplesheet, checkIfExists: true)
           .splitCsv(header: true, sep: ',')
           .filter( ~/.*long.*/ )
           .map { row -> [ row.sampleID, row.SRA_or_FASTQ, row.library_layout ] }
           .set{ samples_list_long_reads }

        Channel.fromPath(params.RNAseq_samplesheet, checkIfExists: true)
           .splitCsv(header: true, sep: ',')
           .filter( ~/.*single.*/ )
           .map { row -> [ row.sampleID, row.SRA_or_FASTQ, row.library_layout ] }
           .set{ samples_list_single_short_reads }

        Channel.fromPath(params.RNAseq_samplesheet, checkIfExists: true)
           .splitCsv(header: true, sep: ',')
           .filter( ~/.*paired.*/ )
           .map { row -> [ row.sampleID, row.SRA_or_FASTQ, row.library_layout ] }
           .set{ samples_list_paired_short_reads }

        // Combined single and paired-end reads
        samples_list_single_short_reads
            .concat(samples_list_paired_short_reads)
            .set { samples_list_short_reads }

        Channel.fromPath(params.protein_samplesheet, checkIfExists: true)
           .splitCsv(header: true, sep: ',')
           .map { row -> [ row.organism, row.filename ] }
           .set{ protein_list }

        evidence_data = generate_evidence_data(
            samples_list_long_reads, 
            samples_list_short_reads, 
            protein_list,
            runtime_flags.run_edta,
            runtime_flags.run_egapx,
            runtime_flags.use_long_reads
        )
    }
    
    // if the workflow parameter is equal to aegis, the outputs generated by the first workflow (generate_evidence_data) must be retrieved manually.
    else if (params.workflow == "aegis") {

      def outdir_1 = params.output_dir
      println "Checking output directory: ${outdir_1}"
      file(outdir_1).mkdirs()

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

      def fake_null_merged_star_stringtie_default_args_unstranded = file("${params.output_dir}/dev_null1")

      if (!fake_null_merged_star_stringtie_default_args_unstranded.exists() || fake_null_merged_star_stringtie_default_args_unstranded.size() == 0) {
          fake_null_merged_star_stringtie_default_args_unstranded.text = "FAKE merged_star_stringtie_default_args_unstranded"
      }

      def fake_null_merged_star_stringtie_unstranded_alt = file("${params.output_dir}/dev_null2")

      if (!fake_null_merged_star_stringtie_unstranded_alt.exists() || fake_null_merged_star_stringtie_unstranded_alt.size() == 0) {
          fake_null_merged_star_stringtie_unstranded_alt.text = "FAKE merged_star_stringtie_unstranded_alt"
      }

      def fake_null_gffcompare_out_star_psiclass_unstranded = file("${params.output_dir}/dev_null3")

      if (!fake_null_gffcompare_out_star_psiclass_unstranded.exists() || fake_null_gffcompare_out_star_psiclass_unstranded.size() == 0) {
          fake_null_gffcompare_out_star_psiclass_unstranded.text = "FAKE gffcompare_out_star_psiclass_unstranded"
      }

      // Load manually the outputs of generate_evidence_data workflow
      if (runtime_flags.run_edta) {
        if (file("${outdir_1}/assembly_masked.EDTA.fasta").exists()) {
          workflow_inputs << tuple("masked_genome.masked_genome", file("${outdir_1}/assembly_masked.EDTA.fasta"))
        } else {
          println "ERROR : File ${outdir_1}/assembly_masked.EDTA.fasta doesn't exists !"
        }
      } else {
        println "EDTA = No, we need this output for Aegis. EXIT."
      }

      if (runtime_flags.use_long_reads) {
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
        workflow_inputs << tuple("merged_star_stringtie.default_args_unstranded", file(fake_null_merged_star_stringtie_default_args_unstranded))
      }
      if (file("${outdir_1}/merged_star_stringtie_unstranded_alt.gtf").exists()) {
          workflow_inputs << tuple("merged_star_stringtie.alt_args_unstranded", file("${outdir_1}/merged_star_stringtie_unstranded_alt.gtf"))
      } else {
        println "WARNING : no RNAseq unstranded data used in this run !"
        workflow_inputs << tuple("merged_star_stringtie.alt_args_unstranded", file(fake_null_merged_star_stringtie_unstranded_alt))
      }
      if (file("${outdir_1}/merged_star_psiclass_unstranded.gtf").exists()) {
          workflow_inputs << tuple("gffcompare_out.star_psiclass_unstranded", file("${outdir_1}/merged_star_psiclass_unstranded.gtf"))
      } else {
        println "WARNING : no RNAseq unstranded data used in this run !"
        workflow_inputs << tuple("gffcompare_out.star_psiclass_unstranded", file(fake_null_gffcompare_out_star_psiclass_unstranded))
      }

      def workflow_inputs_list = workflow_inputs // Keep the list as it is

      println "Files sent to Aegis : ${workflow_inputs_list}"

      aegis(workflow_inputs_list, runtime_flags.run_edta, runtime_flags.use_long_reads) // We send the list directly, not a Channel
  }
}
