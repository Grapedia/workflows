nextflow.enable.dsl = 2

// Include subworkflows
include { generate_evidence_data } from '../subworkflows/generate_evidence_data'
include { aegis } from '../subworkflows/aegis'
include { titan_provenance } from '../modules/titan_provenance'
include { validate_final_annotation } from '../modules/validate_final_annotation'
include { validate_inputs } from '../modules/validate_inputs'

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

def parseCsvLine(String line) {
    def cells = []
    def cell = new StringBuilder()
    boolean quoted = false
    boolean skipNext = false
    line.toList().eachWithIndex { current, i ->
        if (skipNext) {
            skipNext = false
            return
        }
        if (current == '"') {
            if (quoted && i + 1 < line.length() && line.charAt(i + 1) == '"') {
                cell.append('"')
                skipNext = true
            } else {
                quoted = !quoted
            }
        } else if (current == ',' && !quoted) {
            cells << cell.toString().trim()
            cell = new StringBuilder()
        } else {
            cell.append(current)
        }
    }
    cells << cell.toString().trim()
    return cells
}

def samplesheetHasLongReads(samplesheetPath) {
    def rows = file(samplesheetPath).readLines().findAll { line ->
        def trimmed = line.trim()
        trimmed && !trimmed.startsWith('#')
    }

    if (rows.size() < 2) {
        return false
    }

    def header = parseCsvLine(rows[0])
    def layoutIndex = header.indexOf('library_layout')
    if (layoutIndex < 0) {
        error "RNA-seq samplesheet must contain a 'library_layout' column"
    }

    return rows.drop(1).any { row ->
        def columns = parseCsvLine(row)
        columns.size() > layoutIndex && columns[layoutIndex].equalsIgnoreCase('long')
    }
}

def rejectDeprecatedWorkflowParam() {
    if (params.containsKey('workflow')) {
        error "--workflow is no longer supported. TITAN always runs evidence generation followed by Aegis in one graph."
    }
}

def rnaseqLocalFiles(row) {
    def source = row.SRA_or_FASTQ?.toString()?.trim()?.toUpperCase()
    def layout = row.library_layout?.toString()?.trim()?.toLowerCase()
    def sampleID = row.sampleID?.toString()?.trim()

    if (source == 'FASTQ' && layout == 'paired') {
        return [file("${params.RNAseq_data_dir}/${sampleID}_1.fastq.gz"), file("${params.RNAseq_data_dir}/${sampleID}_2.fastq.gz")]
    }
    if (source == 'FASTQ' && layout in ['single', 'long']) {
        return [file("${params.RNAseq_data_dir}/${sampleID}.fastq.gz")]
    }
    if (source == 'FASTA' && layout == 'long') {
        return [file("${params.RNAseq_data_dir}/${sampleID}.fasta")]
    }
    return []
}

workflow TITAN {
    rejectDeprecatedWorkflowParam()
    validateRequiredParams(['output_dir', 'egapx_paramfile', 'RNAseq_samplesheet', 'RNAseq_data_dir', 'protein_samplesheet', 'new_assembly', 'previous_assembly', 'previous_annotations'])
    validateExistingInputFiles(['egapx_paramfile', 'RNAseq_samplesheet', 'protein_samplesheet', 'new_assembly', 'previous_assembly', 'previous_annotations'])
    def has_long_reads = samplesheetHasLongReads(params.RNAseq_samplesheet)
    println "Long-read RNA-seq detected from samplesheet: ${has_long_reads}"

    input_validation = validate_inputs(
        file(params.new_assembly),
        file(params.previous_assembly),
        file(params.previous_annotations),
        file(params.RNAseq_samplesheet),
        params.RNAseq_data_dir,
        file(params.protein_samplesheet),
        file(params.egapx_paramfile),
        params.egapx_paramfile,
        params.egapx_executor,
        params.PSICLASS_vd_option,
        params.PSICLASS_c_option,
        params.egapx_version,
        params.egapx_revision,
        params.container_egapx,
        params.egapx_data_version,
        params.aegis_version,
        params.container_aegis
    )

    new_assembly = input_validation.ok.map { file(params.new_assembly) }
    previous_assembly = input_validation.ok.map { file(params.previous_assembly) }
    previous_annotations = input_validation.ok.map { file(params.previous_annotations) }
    rnaseq_samplesheet = input_validation.ok.map { file(params.RNAseq_samplesheet) }
    protein_samplesheet = input_validation.ok.map { file(params.protein_samplesheet) }
    egapx_paramfile = input_validation.ok.map { file(params.egapx_paramfile) }

    rnaseq_rows = rnaseq_samplesheet
       .splitCsv(header: true, sep: ',')
       .map { row -> [ row.sampleID, row.SRA_or_FASTQ, row.library_layout, rnaseqLocalFiles(row) ] }

    samples_list_long_reads = rnaseq_rows
        .filter { sample_ID, SRA_or_FASTQ, library_layout, local_reads -> library_layout?.toString()?.trim()?.equalsIgnoreCase('long') }

    rnaseq_rows
        .filter { sample_ID, SRA_or_FASTQ, library_layout, local_reads -> library_layout?.toString()?.trim()?.equalsIgnoreCase('single') }
        .concat(
            rnaseq_rows.filter { sample_ID, SRA_or_FASTQ, library_layout, local_reads -> library_layout?.toString()?.trim()?.equalsIgnoreCase('paired') }
        )
        .set { samples_list_short_reads }

    protein_samplesheet
       .splitCsv(header: true, sep: ',')
       .map { row -> [ row.organism, row.filename ] }
       .set{ protein_list }

    evidence_data = generate_evidence_data(
        new_assembly,
        file(params.new_assembly).getName(),
        previous_assembly,
        previous_annotations,
        egapx_paramfile,
        file("${projectDir}/scripts/edta.sh"),
        file("${projectDir}/scripts/Stringtie.sh"),
        file("${projectDir}/scripts/Stringtie_AltCommands.sh"),
        file("${projectDir}/scripts/run_stringtie_transcriptome.sh"),
        file("${projectDir}/scripts/download_sra_fastq.py"),
        file("${projectDir}/scripts/clean_liftoff_gff3_for_agat.py"),
        file("${projectDir}/scripts/clean_protein_fasta_for_BRAKER3.py"),
        file("${projectDir}/scripts/run_braker3_prediction.sh"),
        file("${projectDir}/assets/empty_default.gtf"),
        file("${projectDir}/assets/empty_alt.gtf"),
        file("${projectDir}/assets/empty_psiclass.gtf"),
        params.ena_download_timeout_seconds,
        params.ena_max_download_attempts,
        params.ena_retry_wait_seconds,
        params.ena_verify_md5,
        params.PSICLASS_vd_option,
        params.PSICLASS_c_option,
        params.STAR_genomeSAindexNbases,
        params.STAR_sjdbGTFfile,
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

    validate_final_annotation(
        evidence_data.masked_genome,
        aegis.out.aegis_gff,
        aegis.out.aegis_proteins_all,
        aegis.out.aegis_proteins_main,
        file("${projectDir}/scripts/validate_final_annotation.py")
    )

    titan_provenance(
        has_long_reads,
        params.output_dir,
        params.egapx_version,
        params.egapx_revision,
        params.container_egapx,
        params.egapx_data_version,
        params.aegis_version,
        params.container_aegis,
        workflow.revision ?: '',
        workflow.commitId ?: '',
        workflow.nextflow.version ?: '',
        workflow.profile ?: '',
        workflow.commandLine ?: '',
        new_assembly,
        previous_assembly,
        previous_annotations,
        rnaseq_samplesheet,
        protein_samplesheet,
        egapx_paramfile,
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
        evidence_data.hisat2_stringtie_stranded_default_gtf,
        evidence_data.hisat2_stringtie_stranded_alt_gtf,
        evidence_data.hisat2_stringtie_unstranded_default_gtf,
        evidence_data.hisat2_stringtie_unstranded_alt_gtf,
        evidence_data.long_reads_default_gtf,
        evidence_data.long_reads_alt_gtf,
        aegis.out.aegis_gff,
        aegis.out.aegis_proteins_all,
        aegis.out.aegis_proteins_main,
        evidence_data.edta_versions,
        evidence_data.egapx_versions,
        evidence_data.braker_versions,
        aegis.out.aegis_versions,
        aegis.out.diamond2go_versions,
        validate_final_annotation.out.versions
    )
}
