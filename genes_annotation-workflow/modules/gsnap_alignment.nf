// 2. Aligning RNAseq data on new_assembly_database with GMAP/GSNAP
process gsnap_alignment {

  tag "GSNAP on ${sample_ID}"
  container 'quay.io/biocontainers/gmap:2020.10.14--pl526h2f06484_0'
  containerOptions "--volume $params.outdir/evidence_data/databases/:/databases"
  publishDir "$params.outdir/evidence_data/RNAseq_$stranded_or_unstranded/alignments/new_assembly"
  cpus 4

  input:
    val(new_assembly_database)
    tuple val(sample_ID), val(stranded_or_unstranded), val(paired_or_single), path(reads)

  output:
    tuple val(sample_ID), val(stranded_or_unstranded), val(paired_or_single), file("${sample_ID}_vs_new_assembly.sam")

  script:
    def basename_database = task.ext.prefix ?: "${new_assembly_database.baseName}"
    """
    if [[ $paired_or_single == "paired" ]]
    then
      gsnap --gunzip --nthreads ${task.cpus} --dir /databases --db ${basename_database} --batch 5 --novelsplicing 1 --format sam --output-file ${sample_ID}_vs_new_assembly.sam --nofails ${reads[0]} ${reads[1]}
    elif [[ $paired_or_single == "single" ]]
    then
      gsnap --gunzip --nthreads ${task.cpus} --dir /databases --db ${basename_database} --batch 5 --novelsplicing 1 --format sam --output-file ${sample_ID}_vs_new_assembly.sam --nofails ${reads}
    fi
    """
}
