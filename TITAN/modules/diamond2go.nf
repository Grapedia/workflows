process diamond2go {
  label 'process_aegis'

  tag "Executing diamond2go annotation on $proteins_file_all and $proteins_file_main"
  container params.container_diamond2go
  publishDir "${params.output_dir}/Diamond2GO_outputs", mode: 'copy'
  input:
    path(proteins_file_all)
    path(proteins_file_main)

  output:
    path "final_annotation_proteins_all.diamond2go.tsv", emit: proteins_all_diamond
    path "final_annotation_proteins_main.diamond2go.tsv", emit: proteins_main_diamond
    path "versions.yml", emit: versions



  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")

    echo "[\$DATE] Running diamond2go on $proteins_file_all"
    CMD="perl /Diamond2GO/Diamond2go.pl -d /Diamond2GO/resources/nr_clean_d2go.dmnd -q $proteins_file_all -t protein"
    echo "[\$DATE] Executing: \$CMD"
    perl /Diamond2GO/Diamond2go.pl -d /Diamond2GO/resources/nr_clean_d2go.dmnd -q $proteins_file_all -t protein
    all_diamond=\$(find . -maxdepth 1 -type f -name 'final_annotation_proteins_all-diamond*' | sort | head -n 1)
    test -s "\$all_diamond"
    cp "\$all_diamond" final_annotation_proteins_all.diamond2go.tsv

    echo "[\$DATE] Running diamond2go on $proteins_file_main"
    CMD="perl /Diamond2GO/Diamond2go.pl -d /Diamond2GO/resources/nr_clean_d2go.dmnd -q $proteins_file_main -t protein"
    echo "[\$DATE] Executing: \$CMD"
    perl /Diamond2GO/Diamond2go.pl -d /Diamond2GO/resources/nr_clean_d2go.dmnd -q $proteins_file_main -t protein
    main_diamond=\$(find . -maxdepth 1 -type f -name 'final_annotation_proteins_main-diamond*' | sort | head -n 1)
    test -s "\$main_diamond"
    cp "\$main_diamond" final_annotation_proteins_main.diamond2go.tsv
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """

  stub:
    """
    printf "query\\tsubject\\n" > final_annotation_proteins_all.diamond2go.tsv
    printf "query\\tsubject\\n" > final_annotation_proteins_main.diamond2go.tsv
    printf '"%s":\n  container: "not_recorded"\n' "${task.process}" > versions.yml
    """
}
