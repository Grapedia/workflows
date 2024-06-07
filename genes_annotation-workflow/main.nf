// Path to outdir
params.outdir = "${projectDir}/intermediate_files"

include { prepare_RNAseq_fastq_files } from "./modules/prepare_RNAseq_fastq_files"
include { trimming_fastq } from "./modules/trimming_fastq"
include { gmap_build_database } from "./modules/gmap_build_database"
include { gsnap_alignment } from "./modules/gsnap_alignment"
include { extracting_primary_mapping } from "./modules/extracting_primary_mapping"
include { assembly_transcriptome_stranded } from "./modules/assembly_transcriptome_stranded"
include { assembly_transcriptome_unstranded } from "./modules/assembly_transcriptome_unstranded"
include { conversion_gtf_gff3 } from "./modules/conversion_gtf_gff3"

Channel.fromPath( file(params.RNAseq_samplesheet) )
                    .splitCsv(header: true, sep: ',')
                    .map { row ->
                          def sample_ID = row.sampleID
                          def stranded_or_unstranded = row.stranded_or_unstranded
                          def SRA_or_FASTQ = row.SRA_or_FASTQ
                          def paired_or_single = row.paired_or_single

                          return [sample_ID, stranded_or_unstranded, SRA_or_FASTQ, paired_or_single]
                    }
                    .set{ samples_list }

Channel.fromPath( file(params.protein_samplesheet) )
                    .splitCsv(header: true, sep: ',')
                    .map { row ->
                          def organism = row.organism
                          def filename = row.filename

                          return [organism, filename]
                    }
                    .set{ protein_list }

workflow {
  // ----------------------------------------------------------------------------------------
  //                                Transcriptomes assembly
  // ----------------------------------------------------------------------------------------
  prepare_RNAseq_fastq_files(samples_list)
  trimming_fastq(prepare_RNAseq_fastq_files.out)
  gmap_build_database(params.assemblies_folder,params.new_assembly)
  gsnap_alignment(gmap_build_database.out,trimming_fastq.out)
  extracting_primary_mapping(gsnap_alignment.out) | collect
  extracting_primary_mapping.out.view()
  // assembly_transcriptome_stranded(extracting_primary_mapping.out)
  // assembly_transcriptome_unstranded(extracting_primary_mapping.out)
  // conversion_gtf_gff3(params.assemblies_folder,params.new_assembly, assembly_transcriptome_stranded.out, assembly_transcriptome_unstranded.out)

  // ----------------------------------------------------------------------------------------
  //                                   Protein alignments
  // ----------------------------------------------------------------------------------------
}
