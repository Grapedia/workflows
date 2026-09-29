nextflow.enable.dsl = 2

/*
 * Controlled test of the BRAKER3 inputs given to Mikado.
 *
 * Re-runs only  mikado_prepare -> TransDecoder -> mikado_serialise -> mikado_pick  (the
 * modules and containers of the production pipeline, unchanged) on the evidence tracks
 * already published by a finished TITAN run, then measures the resulting gene set (AGAT
 * statistics, BUSCO on the longest-isoform proteins).  Only the two BRAKER3 tracks differ
 * between the two arms:
 *
 *   --braker_mode raw     augustus.hints.gff3 + genemark.gtf              (production behaviour)
 *   --braker_mode tsebra  braker.gff3         + genemark_supported.gtf    (TSEBRA-filtered)
 *
 * Run it from the project root so that projectDir and nextflow.config match production:
 *   scripts/mikado_input_test/launch_arm.sh raw|tsebra
 * See docs/user/mikado_braker_input_test.md.
 */

include { mikado_prepare; mikado_serialise; mikado_pick } from './modules/mikado'
include { transdecoder_longorfs; transdecoder_predict } from './modules/transdecoder'
include { busco } from './modules/busco'
include { agat_stats } from './modules/agat_stats'

params.titan_output_dir = "${projectDir}/data/titan_prod_out"
params.braker_mode = 'raw'
params.unmasked_genome = "${projectDir}/data/assemblies/T2T_ref.fasta"

process mikado_main_proteins {
  label 'process_low'
  tag "Longest-isoform proteins from the Mikado loci (${params.braker_mode})"
  container params.container_agat
  publishDir "${params.output_dir}/mikado_eval", mode: 'copy'

  input:
    path(genome)
    path(mikado_gff3)

  output:
    path "mikado_proteins_main.fasta", emit: main
    path "mikado_proteins_all.fasta", emit: all

  script:
    """
    set -euo pipefail
    # Mikado loci: gene -> mRNA -> exon/CDS, mRNA IDs are <gene>.<n>.
    # AGAT rejects Mikado's `superlocus` features and needs a genome with fixed-width lines.
    awk -F'\\t' '\$3 != "superlocus"' ${mikado_gff3} > mikado_no_superlocus.gff3
    fold -w 80 ${genome} > genome.wrapped.fa
    agat_sp_extract_sequences.pl -g mikado_no_superlocus.gff3 -f genome.wrapped.fa -p -o mikado_proteins_all.raw.fasta
    # '*' and '-' break HMMER-based tools (BUSCO); keep the longest isoform of every gene.
    # (No python in the AGAT image: plain awk.)
    awk '/^>/ {name = substr(\$1, 2); next} {gsub(/[*-]/, ""); seq[name] = seq[name] \$0}
         END {for (n in seq) print n "\\t" seq[n]}' mikado_proteins_all.raw.fasta > proteins.tsv
    awk -F'\\t' '{print ">" \$1 "\\n" \$2}' proteins.tsv > mikado_proteins_all.fasta
    awk -F'\\t' '{g = \$1; sub(/\\.[0-9]+\$/, "", g); if (length(\$2) > best[g]) {best[g] = length(\$2); id[g] = \$1; sq[g] = \$2}}
         END {for (g in id) print ">" id[g] "\\n" sq[g]}' proteins.tsv > mikado_proteins_main.fasta
    test -s mikado_proteins_main.fasta
    """
}

workflow {
    def te = params.titan_output_dir
    def ev = "${te}/04_evidence"
    def gp = "${ev}/gene_prediction"
    def ta = "${ev}/transcript_assemblies"

    if( !(params.braker_mode in ['raw', 'tsebra']) )
        error "--braker_mode must be 'raw' or 'tsebra', got '${params.braker_mode}'"

    def augustus = file(params.braker_mode == 'tsebra' ? "${gp}/braker.gff3" : "${gp}/augustus.hints.gff3", checkIfExists: true)
    def genemark = file(params.braker_mode == 'tsebra' ? "${gp}/genemark_supported.gtf" : "${gp}/genemark.gtf", checkIfExists: true)
    println "BRAKER3 tracks given to Mikado (${params.braker_mode}): ${augustus.name} + ${genemark.name}"

    // run_hisat2 is false in production: the four HISAT2 slots receive empty placeholders.
    def empty_gtf = file("${projectDir}/data/mikado_input_test/empty.gtf")
    if( !empty_gtf.exists() ) {
        empty_gtf.parent.mkdirs()
        empty_gtf.text = ''
    }

    def genome = file("${ev}/assembly_masked.EDTA.fasta", checkIfExists: true)

    mikado_prepared = mikado_prepare(
        genome,
        augustus,
        genemark,
        file("${ev}/liftoff_previous_annotations.gff3", checkIfExists: true),
        file("${ev}/egapx/egapx.complete.genomic.gff3", checkIfExists: true),
        file("${ta}/merged_star_stringtie_stranded_default.gtf", checkIfExists: true),
        file("${ta}/merged_star_stringtie_stranded_alt.gtf", checkIfExists: true),
        file("${ta}/merged_star_psiclass_stranded.gtf", checkIfExists: true),
        file("${ta}/merged_star_psiclass_unstranded.gtf", checkIfExists: true),
        file("${ta}/merged_star_stringtie_unstranded_default.gtf", checkIfExists: true),
        file("${ta}/merged_star_stringtie_unstranded_alt.gtf", checkIfExists: true),
        empty_gtf, empty_gtf, empty_gtf, empty_gtf,
        file("${ta}/merged_minimap2_stringtie_long_reads_default.gtf", checkIfExists: true),
        file("${ta}/merged_minimap2_stringtie_long_reads_alt.gtf", checkIfExists: true),
        file("${te}/03_additional_annotations/flair/flair_isoforms.gtf", checkIfExists: true),
        file("${te}/03_additional_annotations/helixer/helixer.gff3", checkIfExists: true),
        file("${projectDir}/scripts/make_mikado_list.py")
    )

    orfs_dir = transdecoder_longorfs(mikado_prepared.fasta)
    orfs = transdecoder_predict(mikado_prepared.fasta, orfs_dir.longorfs_dir)
    db = mikado_serialise(
        mikado_prepared.config,
        mikado_prepared.fasta,
        orfs.bed,
        file("${projectDir}/assets/mikado_pandas_sqlalchemy_sitecustomize.py")
    )
    picked = mikado_pick(genome, mikado_prepared.config, mikado_prepared.gtf, db.database)

    proteins = mikado_main_proteins(file(params.unmasked_genome, checkIfExists: true), picked.gff3)
    busco(proteins.main)
    agat_stats(picked.gff3)
}
