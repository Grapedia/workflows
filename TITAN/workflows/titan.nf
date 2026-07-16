nextflow.enable.dsl = 2

// Include subworkflows
include { generate_evidence_data } from '../subworkflows/generate_evidence_data'
include { aegis } from '../subworkflows/aegis'

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
    def allowedWorkflows = ['generate_evidence_data', 'aegis', 'all']
    if (!allowedWorkflows.contains(params.workflow)) {
        error "Invalid --workflow '${params.workflow}'. Allowed values: ${allowedWorkflows.join(', ')}"
    }
}

def samplesheetHasLongReads(samplesheetPath) {
    def rows = file(samplesheetPath).readLines().findAll { line ->
        def trimmed = line.trim()
        trimmed && !trimmed.startsWith('#')
    }

    if (rows.size() < 2) {
        return false
    }

    def header = rows[0].split(',', -1).collect { it.trim() }
    def layoutIndex = header.indexOf('library_layout')
    if (layoutIndex < 0) {
        error "RNA-seq samplesheet must contain a 'library_layout' column"
    }

    return rows.drop(1).any { row ->
        def columns = row.split(',', -1).collect { it.trim() }
        columns.size() > layoutIndex && columns[layoutIndex].equalsIgnoreCase('long')
    }
}

workflow TITAN {
    validateWorkflowName()
    validateRequiredParams(['output_dir', 'egapx_paramfile', 'RNAseq_samplesheet', 'protein_samplesheet', 'new_assembly', 'previous_assembly', 'previous_annotations'])
    validateExistingInputFiles(['egapx_paramfile', 'RNAseq_samplesheet', 'protein_samplesheet', 'new_assembly', 'previous_assembly', 'previous_annotations'])
    def has_long_reads = samplesheetHasLongReads(params.RNAseq_samplesheet)
    println "Long-read RNA-seq detected from samplesheet: ${has_long_reads}"

    if (params.workflow == "generate_evidence_data" || params.workflow == "all") {
        Channel.fromPath(params.RNAseq_samplesheet, checkIfExists: true)
           .splitCsv(header: true, sep: ',')
           .filter { row -> row.library_layout?.toString()?.trim()?.equalsIgnoreCase('long') }
           .map { row -> [ row.sampleID, row.SRA_or_FASTQ, row.library_layout ] }
           .set{ samples_list_long_reads }

        Channel.fromPath(params.RNAseq_samplesheet, checkIfExists: true)
           .splitCsv(header: true, sep: ',')
           .filter { row -> row.library_layout?.toString()?.trim()?.equalsIgnoreCase('single') }
           .map { row -> [ row.sampleID, row.SRA_or_FASTQ, row.library_layout ] }
           .set{ samples_list_single_short_reads }

        Channel.fromPath(params.RNAseq_samplesheet, checkIfExists: true)
           .splitCsv(header: true, sep: ',')
           .filter { row -> row.library_layout?.toString()?.trim()?.equalsIgnoreCase('paired') }
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
            has_long_reads
        )

        if (params.workflow == "all") {
            println "Running Aegis from generated named evidence; EDTA masked genome is passed as a direct channel input."

            aegis(
                evidence_data.masked_genome,
                evidence_data.braker_augustus_gff3,
                evidence_data.braker_genemark_gtf,
                evidence_data.liftoff_gff3,
                evidence_data.star_stringtie_stranded_default_gtf,
                evidence_data.star_stringtie_stranded_alt_gtf,
                evidence_data.star_psiclass_stranded_gtf,
                evidence_data.star_psiclass_unstranded_gtf,
                evidence_data.star_stringtie_unstranded_default_gtf,
                evidence_data.star_stringtie_unstranded_alt_gtf,
                evidence_data.long_reads_default_gtf,
                evidence_data.long_reads_alt_gtf,
                has_long_reads
            )
        }
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

      def missing_required_inputs = []
      def requiredEvidenceFile = { filename ->
        def candidate = file("${outdir_1}/${filename}")
        if (candidate.exists()) {
          return candidate
        }
        missing_required_inputs << filename
        return candidate
      }
      def optionalEvidenceFile = { filename, fallback ->
        def candidate = file("${outdir_1}/${filename}")
        if (candidate.exists()) {
          return candidate
        }
        println "WARNING : optional evidence ${filename} not found in ${outdir_1}; using ${fallback.name}"
        return fallback
      }

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

      def fake_null_long_reads_default = file("${params.output_dir}/dev_null_long1")

      if (!fake_null_long_reads_default.exists() || fake_null_long_reads_default.size() == 0) {
          fake_null_long_reads_default.text = "FAKE merged_long_reads_default"
      }

      def fake_null_long_reads_alt = file("${params.output_dir}/dev_null_long2")

      if (!fake_null_long_reads_alt.exists() || fake_null_long_reads_alt.size() == 0) {
          fake_null_long_reads_alt.text = "FAKE merged_long_reads_alt"
      }

      // Load manually the outputs of generate_evidence_data workflow
      def masked_genome = requiredEvidenceFile.call("assembly_masked.EDTA.fasta")

      def long_reads_default = fake_null_long_reads_default
      def long_reads_alt = fake_null_long_reads_alt
      if (has_long_reads) {
        long_reads_default = requiredEvidenceFile.call("merged_minimap2_stringtie_long_reads_default.gtf")
        long_reads_alt = requiredEvidenceFile.call("merged_minimap2_stringtie_long_reads_alt.gtf")
      }

      def braker_augustus_gff = requiredEvidenceFile.call("augustus.hints.gff3")
      def braker_genemark_gtf = requiredEvidenceFile.call("genemark.gtf")
      def liftoff_annotation = requiredEvidenceFile.call("liftoff_previous_annotations.gff3")
      def star_stringtie_default_stranded = requiredEvidenceFile.call("merged_star_stringtie_stranded_default.gtf")
      def star_stringtie_alt_stranded = requiredEvidenceFile.call("merged_star_stringtie_stranded_alt.gtf")
      def star_psiclass_stranded = requiredEvidenceFile.call("merged_star_psiclass_stranded.gtf")
      def star_stringtie_default_unstranded = optionalEvidenceFile.call("merged_star_stringtie_unstranded_default.gtf", fake_null_merged_star_stringtie_default_args_unstranded)
      def star_stringtie_alt_unstranded = optionalEvidenceFile.call("merged_star_stringtie_unstranded_alt.gtf", fake_null_merged_star_stringtie_unstranded_alt)
      def star_psiclass_unstranded = optionalEvidenceFile.call("merged_star_psiclass_unstranded.gtf", fake_null_gffcompare_out_star_psiclass_unstranded)

      if (missing_required_inputs) {
        error "Missing required Aegis evidence file(s) in ${outdir_1}:\n  ${missing_required_inputs.join('\n  ')}"
      }

      println "Named evidence sent to Aegis: masked_genome=${masked_genome}, braker_augustus_gff=${braker_augustus_gff}, braker_genemark_gtf=${braker_genemark_gtf}, liftoff_annotation=${liftoff_annotation}"

      aegis(
          masked_genome,
          braker_augustus_gff,
          braker_genemark_gtf,
          liftoff_annotation,
          star_stringtie_default_stranded,
          star_stringtie_alt_stranded,
          star_psiclass_stranded,
          star_psiclass_unstranded,
          star_stringtie_default_unstranded,
          star_stringtie_alt_unstranded,
          long_reads_default,
          long_reads_alt,
          has_long_reads
      )
  }
}
