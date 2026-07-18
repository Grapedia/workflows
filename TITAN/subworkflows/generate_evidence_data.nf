nextflow.enable.dsl = 2

// Include various processing modules
include { prepare_RNAseq_fastq_files_short } from "../modules/prepare_RNAseq_fastq_files_short"
include { prepare_RNAseq_fastq_files_long } from "../modules/prepare_RNAseq_fastq_files_long"
include { trimming_fastq } from "../modules/trimming_fastq"
include { liftoff_annotations } from "../modules/liftoff_annotations"
include { egapx } from "../modules/egapx"
include { clean_liftoff_gff3_for_agat } from "../modules/clean_liftoff_gff3_for_agat"
include { agat_convert_gff3_to_cds_fasta } from "../modules/agat_convert_gff3_to_cds_fasta"
include { salmon_index } from "../modules/salmon_index"
include { salmon_strand_inference } from "../modules/salmon_strand_inference"
include { star_genome_indices } from "../modules/star_genome_indices"
include { star_alignment } from "../modules/star_alignment"
include { hisat2_genome_indices } from "../modules/hisat2_genome_indices"
include { hisat2_alignment } from "../modules/hisat2_alignment"
include { minimap2_genome_indices } from "../modules/minimap2_genome_indices"
include { minimap2_alignment } from "../modules/minimap2_alignment"
include { assembly_transcriptome_star_psiclass } from "../modules/assembly_transcriptome_star_psiclass"
include { assembly_transcriptome_star_stringtie } from "../modules/assembly_transcriptome_star_stringtie"
include { assembly_transcriptome_hisat2_stringtie } from "../modules/assembly_transcriptome_hisat2_stringtie"
include { gffcompare } from "../modules/gffcompare"
include { assembly_transcriptome_minimap2_stringtie } from "../modules/assembly_transcriptome_minimap2_stringtie"
include { Stringtie_merging_short_reads_STAR } from "../modules/Stringtie_merging_short_reads_STAR"
include { Stringtie_merging_short_reads_hisat2 } from "../modules/Stringtie_merging_short_reads_hisat2"
include { Stringtie_merging_long_reads } from "../modules/Stringtie_merging_long_reads"
include { EDTA } from "../modules/EDTA"
include { braker3_prediction } from "../modules/braker3_prediction"
include { braker3_prediction_with_long_reads } from "../modules/braker3_prediction_with_long_reads"
include { normalize_protein_fastas } from "../modules/normalize_protein_fastas"
include { empty_long_read_evidence } from "../modules/empty_long_read_evidence"
include { flair_isoforms; flair_merge_isoforms; flair_empty_evidence } from "../modules/flair"

def collectTranscriptGtfs(transcript_channel, String target_strand, Integer file_index, fallback_gtf = null) {
    def collected = transcript_channel
      .filter { row -> target_strand == 'stranded' ? row[3] != 'unstranded' : row[3] == 'unstranded' }
      .map { row -> row[file_index] }
      .collect()

    return fallback_gtf == null ? collected : collected.ifEmpty(fallback_gtf)
}

def collectPsiclassGtfs(psiclass_channel, String target_strand, fallback_gtf = null) {
    def collected = psiclass_channel
      .filter { sample_ID, gtf_file, strand_type -> target_strand == 'stranded' ? strand_type != 'unstranded' : strand_type == 'unstranded' }
      .map { sample_ID, gtf_file, strand_type -> gtf_file }
      .collect()

    return fallback_gtf == null ? collected : collected.ifEmpty(fallback_gtf)
}

workflow generate_evidence_data {
    take:
        new_assembly
        new_assembly_name
        previous_assembly
        previous_annotations_file
        egapx_paramfile
        edta_script
        stringtie_script
        stringtie_alt_script
        stringtie_transcriptome_script
        download_sra_fastq_script
        clean_liftoff_gff3_script
        clean_protein_script
        braker3_runner_script
        empty_default_gtf
        empty_alt_gtf
        empty_psiclass_gtf
        ena_download_timeout_seconds
        ena_max_download_attempts
        ena_retry_wait_seconds
        ena_verify_md5
        psiclass_vd_option
        psiclass_c_option
        star_genome_sa_index_nbases
        star_sjdb_gtf_file
        // RNA-seq tuples:
        //   long:  tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path(local_reads)
        //   short: tuple val(sample_ID), val(SRA_or_FASTQ), val(library_layout), path(local_reads)
        // Protein tuples: tuple val(organism), val(filename)
        samples_list_long_reads
        samples_list_short_reads
        protein_list
        has_long_reads

    main:
        // Prepare RNAseq short reads for processing
        short_reads_prepared = prepare_RNAseq_fastq_files_short(
            samples_list_short_reads,
            download_sra_fastq_script,
            ena_download_timeout_seconds,
            ena_max_download_attempts,
            ena_retry_wait_seconds,
            ena_verify_md5
        )

        // Trim Illumina short reads
        trimmed_reads = trimming_fastq(short_reads_prepared.prepared_fastqs)

        // Lift over previous annotations to new assembly
        previous_annotations = liftoff_annotations(
            new_assembly,
            previous_assembly,
            previous_annotations_file
        )

        // EGAPx annotation pipeline on new assembly
        egapx_annotations = egapx(egapx_paramfile)

        cleaned_liftoff_gff3 = clean_liftoff_gff3_for_agat(
            previous_annotations.liftoff_previous_annotations,
            clean_liftoff_gff3_script
        )

        // Convert cleaned GFF3 to CDS FASTA for Salmon strand inference
        gff_cds = agat_convert_gff3_to_cds_fasta(
            new_assembly,
            cleaned_liftoff_gff3.cleaned_gff3
        )

        // Salmon index and strand inference.
        // Emits tuple val(sample_ID), val(library_layout), path(read_1), path(read_2), val(strand_type).
        salmon_index_result = salmon_index(gff_cds.cds_fasta)
        strand_inference = salmon_strand_inference(trimmed_reads.trimmed_reads, salmon_index_result.index)

        // Salmon emits the inferred strand as a value and keeps logs/classification as debug outputs.
        salmon_output_processed = strand_inference.strand_inference_result.map { sample_ID, library_layout, read_1, read_2, strand_info ->
            return [sample_ID, library_layout, read_1, read_2, strand_info]
        }

        // Process STAR alignments.
        // Emits tuple val(sample_ID), path(bam_file), val(strand_type).
        star_indices = star_genome_indices(
            new_assembly,
            new_assembly_name,
            star_genome_sa_index_nbases,
            star_sjdb_gtf_file
        )
        star_aligned = star_alignment(star_indices.index, salmon_output_processed)

        star_aligned.samples_aligned
          .map { sample_ID, bam_file, strand_type -> bam_file }
          .collect()
          .set{ concat_star_bams_BRAKER3 }

        // Process HISAT2 alignments
        hisat2_indices = hisat2_genome_indices(
            new_assembly,
            new_assembly_name
        )
        hisat2_aligned = hisat2_alignment(hisat2_indices.index, salmon_output_processed, new_assembly_name)

        // Process hisat2 assemblies
        hisat2_assemblies_stringtie = assembly_transcriptome_hisat2_stringtie(
            hisat2_aligned.samples_aligned,
            stringtie_script,
            stringtie_alt_script,
            stringtie_transcriptome_script
        )

        hisat2_stranded_default_gtfs = collectTranscriptGtfs(hisat2_assemblies_stringtie.hisat2_stringtie_transcriptomes, 'stranded', 1)
        hisat2_stranded_alt_gtfs = collectTranscriptGtfs(hisat2_assemblies_stringtie.hisat2_stringtie_transcriptomes, 'stranded', 2)
        hisat2_unstranded_default_gtfs = collectTranscriptGtfs(hisat2_assemblies_stringtie.hisat2_stringtie_transcriptomes, 'unstranded', 1, empty_default_gtf)
        hisat2_unstranded_alt_gtfs = collectTranscriptGtfs(hisat2_assemblies_stringtie.hisat2_stringtie_transcriptomes, 'unstranded', 2, empty_alt_gtf)

        merged_hisat2_stringtie = Stringtie_merging_short_reads_hisat2(
            hisat2_stranded_default_gtfs,
            hisat2_stranded_alt_gtfs,
            hisat2_unstranded_default_gtfs,
            hisat2_unstranded_alt_gtfs
        )

        if (!has_long_reads) {
            empty_long_reads = empty_long_read_evidence(empty_default_gtf, empty_alt_gtf)
            merged_long_reads_default_args_gtf = empty_long_reads.default_args_gtf
            merged_long_reads_alt_args_gtf = empty_long_reads.alt_args_gtf
            empty_flair = flair_empty_evidence(empty_default_gtf)
            merged_flair_isoforms_gtf = empty_flair.gtf
            merged_flair_isoforms_fasta = empty_flair.fasta
        }

        // Process minimap2 alignments for long reads
        if (has_long_reads) {
            // Prepare long reads (if any) for processing
            long_reads_prepared = prepare_RNAseq_fastq_files_long(
                samples_list_long_reads,
                download_sra_fastq_script,
                ena_download_timeout_seconds,
                ena_max_download_attempts,
                ena_retry_wait_seconds,
                ena_verify_md5
            )

            minimap2_indices = minimap2_genome_indices(
                new_assembly,
                new_assembly_name
            )
            minimap2_aligned = minimap2_alignment(minimap2_indices.index, long_reads_prepared.prepared_fastqs)

            minimap2_aligned.samples_aligned
              .map { sample_ID, bam_file -> bam_file }
              .collect()
              .set { concat_minimap2_bams_BRAKER3 }

            // Process long read assemblies
            long_reads_assemblies = assembly_transcriptome_minimap2_stringtie(
                minimap2_aligned.samples_aligned,
                stringtie_script,
                stringtie_alt_script,
                stringtie_transcriptome_script
            )

            long_reads_assemblies.minimap2_stringtie_transcriptomes
              .map { sample_ID, default_gtf, alt_gtf -> default_gtf }
              .collect()
              .set { minimap2_default_gtfs }

            long_reads_assemblies.minimap2_stringtie_transcriptomes
              .map { sample_ID, default_gtf, alt_gtf -> alt_gtf }
              .collect()
              .set { minimap2_alt_gtfs }

            merged_long_reads = Stringtie_merging_long_reads(minimap2_default_gtfs, minimap2_alt_gtfs)
            merged_long_reads_default_args_gtf = merged_long_reads.default_args_gtf
            merged_long_reads_alt_args_gtf = merged_long_reads.alt_args_gtf

            flair_results = flair_isoforms(
                long_reads_prepared.prepared_fastqs,
                new_assembly,
                previous_annotations.liftoff_previous_annotations
            )
            flair_merged = flair_merge_isoforms(
                flair_results.isoforms.map { sample_ID, gtf, fasta -> gtf }.collect(),
                flair_results.isoforms.map { sample_ID, gtf, fasta -> fasta }.collect()
            )
            merged_flair_isoforms_gtf = flair_merged.gtf
            merged_flair_isoforms_fasta = flair_merged.fasta
        }

        // Process short read assemblies
        star_assemblies_stringtie = assembly_transcriptome_star_stringtie(
            star_aligned.samples_aligned,
            stringtie_script,
            stringtie_alt_script,
            stringtie_transcriptome_script
        )
        star_assemblies_psiclass = assembly_transcriptome_star_psiclass(
            star_aligned.samples_aligned,
            psiclass_vd_option,
            psiclass_c_option
        )

        star_stranded_default_gtfs = collectTranscriptGtfs(star_assemblies_stringtie.star_stringtie_transcriptomes, 'stranded', 1)
        star_stranded_alt_gtfs = collectTranscriptGtfs(star_assemblies_stringtie.star_stringtie_transcriptomes, 'stranded', 2)
        star_unstranded_default_gtfs = collectTranscriptGtfs(star_assemblies_stringtie.star_stringtie_transcriptomes, 'unstranded', 1, empty_default_gtf)
        star_unstranded_alt_gtfs = collectTranscriptGtfs(star_assemblies_stringtie.star_stringtie_transcriptomes, 'unstranded', 2, empty_alt_gtf)

        // Merge assemblies
        merged_star_stringtie = Stringtie_merging_short_reads_STAR(
            star_stranded_default_gtfs,
            star_stranded_alt_gtfs,
            star_unstranded_default_gtfs,
            star_unstranded_alt_gtfs
        )

        star_psiclass_stranded_gtfs = collectPsiclassGtfs(star_assemblies_psiclass.psiclass_assemblies, 'stranded')
        star_psiclass_unstranded_gtfs = collectPsiclassGtfs(star_assemblies_psiclass.psiclass_assemblies, 'unstranded', empty_psiclass_gtf)

        // GFFcompare to merge PsiCLASS transcriptomes
        gffcompare_out = gffcompare(star_psiclass_stranded_gtfs, star_psiclass_unstranded_gtfs)

        // EDTA is mandatory: Aegis requires the hard-masked genome.
        edta_results = EDTA(
            new_assembly,
            new_assembly_name,
            edta_script
        )

        normalized_proteins = normalize_protein_fastas(
            protein_list.map { organism, filename -> [organism, file(filename)] },
            clean_protein_script
        )
        protein_fastas = normalized_proteins.normalized_fastas
          .map { organism, protein_fasta -> protein_fasta }
          .collect()
          .ifEmpty { error "protein_samplesheet must contain at least one protein FASTA" }

        // Run BRAKER3
        if (has_long_reads) {
            braker3_results = braker3_prediction_with_long_reads(
                new_assembly,
                protein_fastas,
                concat_star_bams_BRAKER3,
                concat_minimap2_bams_BRAKER3,
                braker3_runner_script
            )
        } else {
            braker3_results = braker3_prediction(
                new_assembly,
                protein_fastas,
                concat_star_bams_BRAKER3,
                braker3_runner_script
            )
        }

    emit:
        masked_genome = edta_results.masked_genome
        egapx_gff3 = egapx_annotations.gff3
        egapx_gtf = egapx_annotations.gtf
        egapx_proteins = egapx_annotations.proteins
        egapx_cds = egapx_annotations.cds
        egapx_transcripts = egapx_annotations.transcripts
        egapx_annotated_genome_asn = egapx_annotations.annotated_genome_asn
        egapx_output_dir = egapx_annotations.output_dir
        egapx_versions = egapx_annotations.versions
        edta_versions = edta_results.versions
        liftoff_gff3 = previous_annotations.liftoff_previous_annotations
        liftoff_unmapped_features = previous_annotations.unmapped_features
        braker_augustus_gff3 = braker3_results.augustus_gff
        braker_genemark_gtf = braker3_results.genemark_gtf
        braker_genemark_supported_gtf = braker3_results.genemark_supported_gtf
        braker_gff3 = braker3_results.braker_gff
        braker_versions = braker3_results.versions
        star_stringtie_stranded_default_gtf = merged_star_stringtie.default_args_stranded
        star_stringtie_stranded_alt_gtf = merged_star_stringtie.alt_args_stranded
        star_stringtie_unstranded_default_gtf = merged_star_stringtie.default_args_unstranded
        star_stringtie_unstranded_alt_gtf = merged_star_stringtie.alt_args_unstranded
        hisat2_stringtie_stranded_default_gtf = merged_hisat2_stringtie.default_args_stranded
        hisat2_stringtie_stranded_alt_gtf = merged_hisat2_stringtie.alt_args_stranded
        hisat2_stringtie_unstranded_default_gtf = merged_hisat2_stringtie.default_args_unstranded
        hisat2_stringtie_unstranded_alt_gtf = merged_hisat2_stringtie.alt_args_unstranded
        star_psiclass_stranded_gtf = gffcompare_out.star_psiclass_stranded
        star_psiclass_unstranded_gtf = gffcompare_out.star_psiclass_unstranded
        long_reads_default_gtf = merged_long_reads_default_args_gtf
        long_reads_alt_gtf = merged_long_reads_alt_args_gtf
        flair_isoforms_gtf = merged_flair_isoforms_gtf
        flair_isoforms_fasta = merged_flair_isoforms_fasta
        fastp_json_reports = trimmed_reads.fastp_json.collect()
}
