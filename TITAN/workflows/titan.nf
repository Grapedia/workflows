nextflow.enable.dsl = 2

// Include subworkflows
include { generate_evidence_data } from '../subworkflows/generate_evidence_data'
include { aegis } from '../subworkflows/aegis'
include { titan_provenance } from '../modules/titan_provenance'

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

def rejectDeprecatedWorkflowParam() {
    if (params.containsKey('workflow')) {
        error "--workflow is no longer supported. TITAN always runs evidence generation followed by Aegis in one graph."
    }
}

workflow TITAN {
    rejectDeprecatedWorkflowParam()
    validateRequiredParams(['output_dir', 'egapx_paramfile', 'RNAseq_samplesheet', 'protein_samplesheet', 'new_assembly', 'previous_assembly', 'previous_annotations'])
    validateExistingInputFiles(['egapx_paramfile', 'RNAseq_samplesheet', 'protein_samplesheet', 'new_assembly', 'previous_assembly', 'previous_annotations'])
    def has_long_reads = samplesheetHasLongReads(params.RNAseq_samplesheet)
    println "Long-read RNA-seq detected from samplesheet: ${has_long_reads}"

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

    println "Running Aegis from generated named evidence; EDTA masked genome is passed as a direct channel input."

    aegis(
        evidence_data.masked_genome,
        evidence_data.braker_augustus_gff3,
        evidence_data.braker_genemark_gtf,
        evidence_data.liftoff_gff3,
        evidence_data.egapx_gff3,
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

    titan_provenance(
        has_long_reads,
        file(params.new_assembly),
        file(params.previous_assembly),
        file(params.previous_annotations),
        file(params.RNAseq_samplesheet),
        file(params.protein_samplesheet),
        file(params.egapx_paramfile),
        evidence_data.masked_genome,
        evidence_data.liftoff_gff3,
        evidence_data.egapx_gff3,
        evidence_data.braker_augustus_gff3,
        evidence_data.braker_genemark_gtf,
        evidence_data.star_stringtie_stranded_default_gtf,
        evidence_data.star_stringtie_stranded_alt_gtf,
        evidence_data.star_psiclass_stranded_gtf,
        evidence_data.star_psiclass_unstranded_gtf,
        evidence_data.star_stringtie_unstranded_default_gtf,
        evidence_data.star_stringtie_unstranded_alt_gtf,
        evidence_data.long_reads_default_gtf,
        evidence_data.long_reads_alt_gtf,
        aegis.out.aegis_gff,
        aegis.out.aegis_proteins_all,
        aegis.out.aegis_proteins_main
    )
}
