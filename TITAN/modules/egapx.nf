process egapx {
  label 'process_prediction'

  tag "Executing NCBI egapx gene annotation pipeline ..."
  // Intentionally no process container: this wrapper launches the official EGAPx
  // nested Nextflow runner, which must access the host Docker/Apptainer runtime.
  publishDir "${params.output_dir}/04_evidence/egapx", mode: 'copy'
  input:
    path egapx_paramfile

  output:
    path "egapx.complete.genomic.gff3", emit: gff3
    path "egapx.complete.genomic.gtf", emit: gtf
    path "egapx.complete.proteins.faa", emit: proteins
    path "egapx.complete.cds.fna", emit: cds
    path "egapx.complete.transcripts.fna", emit: transcripts
    path "egapx.annotated_genome.asn", emit: annotated_genome_asn
    path "egapx_out", emit: output_dir
    path "versions.yml", emit: versions

  script:
    """
    set -euo pipefail
    DATE=\$(date "+%Y-%m-%d %H:%M:%S")
    TASK_PWD="\$PWD"

    EGAPX_VERSION="${params.egapx_version}"
    EGAPX_REVISION="${params.egapx_revision}"
    EGAPX_CONTAINER="${params.container_egapx}"
    EGAPX_EXECUTOR="${params.egapx_executor}"
    EGAPX_PYTHON="${params.egapx_python}"
    EGAPX_DATA_VERSION="${params.egapx_data_version}"
    EGAPX_RUNNER_DIR="${params.egapx_runner_dir}"
    EGAPX_LOCAL_CACHE_DIR="${params.egapx_local_cache_dir}"
    EGAPX_CONFIG_DIR="${params.egapx_config_dir}"
    EGAPX_WORK_DIR="${params.egapx_work_dir}"

    for required_command in "\$EGAPX_PYTHON" curl tar "\$EGAPX_EXECUTOR"; do
      if ! command -v "\$required_command" >/dev/null 2>&1; then
        echo "ERROR: EGAPx wrapper requires '\$required_command' on the host PATH because it launches a nested Nextflow workflow." >&2
        exit 1
      fi
    done

    # Resolve the paramfile to an absolute path now, while we're still in the
    # task's own work dir -- EGAPX_WORK_DIR mode below `cd`s elsewhere before
    # invoking egapx.py, so a relative/symlinked path wouldn't resolve there.
    EGAPX_PARAMFILE="\$(readlink -f "${egapx_paramfile}")"

    # When EGAPX_WORK_DIR is set, EGAPx's own work dir/output dir/.nextflow
    # session live at a fixed path shared across every attempt of this
    # process, instead of inside this (ephemeral, per-attempt) task work dir.
    # That's what makes resuming a partially-completed nested run possible:
    # Nextflow's own resume is keyed on the .nextflow session found in the
    # directory the nested `nextflow run` is launched from, and on the -w
    # work dir passed to it -- both must be stable across attempts, which the
    # TITAN task work dir (a new hash whenever egapx_paramfile's content
    # changes, e.g. after editing out a bad sample) is not.
    if [[ -n "\$EGAPX_WORK_DIR" && "\$EGAPX_WORK_DIR" != "false" ]]; then
      PERSIST_DIR="\$EGAPX_WORK_DIR"
    else
      PERSIST_DIR="\$TASK_PWD"
    fi
    mkdir -p "\$PERSIST_DIR/egapx_runner" "\$PERSIST_DIR/egapx_work" "\$PERSIST_DIR/egapx_out" "\$PERSIST_DIR/egapx_tmp"

    # A minimap2_wnode sub-task once failed with ENOENT from 'sort' writing to
    # \$TMPDIR (the shared \${projectDir}/.tmp from slurm_apptainer.config)
    # inside its nested singularity container. EGAPx's nested Nextflow session
    # binds only its own task workDir into that container, so an external
    # TMPDIR isn't guaranteed visible there -- though re-testing with the
    # actual production apptainer engine (1.4.0-rc.2) didn't reproduce the gap,
    # so the true trigger (possibly a transient host/storage hiccup instead)
    # is unconfirmed. Redirecting TMPDIR under PERSIST_DIR costs nothing and
    # keeps it inside the tree EGAPx already binds into every nested container.
    export TMPDIR="\$PERSIST_DIR/egapx_tmp"

    if [[ -e "\$PERSIST_DIR/egapx_runner/ui/egapx.py" ]]; then
      echo "[\$DATE] Reusing EGAPx runner already staged at \$PERSIST_DIR/egapx_runner"
    elif [[ -n "\$EGAPX_RUNNER_DIR" && "\$EGAPX_RUNNER_DIR" != "false" ]]; then
      echo "[\$DATE] Using pre-staged EGAPx runner from \$EGAPX_RUNNER_DIR"
      cp -R "\$EGAPX_RUNNER_DIR"/. "\$PERSIST_DIR/egapx_runner/"
    else
      echo "[\$DATE] Downloading EGAPx runner revision \$EGAPX_REVISION. For strict offline reproducibility, provide --egapx_runner_dir."
      curl -fsSL "https://github.com/ncbi/egapx/archive/refs/tags/\${EGAPX_REVISION}.tar.gz" \\
        | tar -xz --strip-components=1 -C "\$PERSIST_DIR/egapx_runner"
    fi

    if [[ -n "\$EGAPX_LOCAL_CACHE_DIR" && "\$EGAPX_LOCAL_CACHE_DIR" != "false" ]]; then
      mkdir -p "\$EGAPX_LOCAL_CACHE_DIR"
      local_cache="\$EGAPX_LOCAL_CACHE_DIR"
    else
      mkdir -p "\$PERSIST_DIR/egapx_local_cache"
      local_cache="\$PERSIST_DIR/egapx_local_cache"
    fi

    printf "process.container = '%s'\\n" "\$EGAPX_CONTAINER" > "\$PERSIST_DIR/egapx_runner/ui/assets/config/docker_image.config"

    CMD=(
      "\$EGAPX_PYTHON" "\$PERSIST_DIR/egapx_runner/ui/egapx.py"
      "\$EGAPX_PARAMFILE"
      -e "\$EGAPX_EXECUTOR"
      -w "\$PERSIST_DIR/egapx_work"
      -o "\$PERSIST_DIR/egapx_out"
      -lc "\$local_cache"
      -dv "\$EGAPX_DATA_VERSION"
    )
    if [[ -n "\$EGAPX_CONFIG_DIR" && "\$EGAPX_CONFIG_DIR" != "false" ]]; then
      if [[ ! -d "\$EGAPX_CONFIG_DIR" ]]; then
        echo "ERROR: EGAPx config directory does not exist: \$EGAPX_CONFIG_DIR" >&2
        exit 1
      fi
      CMD+=(-c "\$EGAPX_CONFIG_DIR")
    fi

    # egapx.py has no --resume flag of its own: on failure it only *writes*
    # an `egapx_out/nextflow/resume.sh` helper for a human to run by hand. To
    # resume automatically we replicate what that helper does: run egapx.py
    # once with -n (dry-run) to refresh egapx_out/nextflow/run_params.yaml
    # from the current paramfile (picking up e.g. samples removed/added since
    # the failed attempt) without touching the nested work dir, then execute
    # the nf_cmd it prints ourselves with -resume appended. Both the dry-run
    # and the real (re)invocation are launched from PERSIST_DIR so the nested
    # `.nextflow` session/cache directory -- which Nextflow looks for relative
    # to the launch directory, not the -w work dir -- is found.
    if [[ -d "\$PERSIST_DIR/.nextflow" ]]; then
      echo "[\$DATE] Found a previous EGAPx nested-run session in \$PERSIST_DIR - resuming instead of restarting from scratch."
      # -n/--dry-run prints the underlying `nextflow run ...` command instead
      # of executing it, but egapx.py still unconditionally prints other
      # diagnostics before it (e.g. "Checking for BUSCO lineage...") and, at
      # default verbosity, a reads summary *after* it with no trailing
      # newline. -q silences that trailing summary so `tail -n 1` reliably
      # isolates just the nf_cmd line.
      NF_CMD="\$(cd "\$PERSIST_DIR" && "\${CMD[@]}" -n -q | tail -n 1)"
      echo "[\$DATE] Resuming: \$NF_CMD -resume"
      ( cd "\$PERSIST_DIR" && NXF_WORK="\$PERSIST_DIR/egapx_work" eval "\$NF_CMD -resume" )
    else
      echo "[\$DATE] Executing: \${CMD[*]}"
      ( cd "\$PERSIST_DIR" && "\${CMD[@]}" )
    fi

    copy_required_egapx_output() {
      local source_path="\$1"
      local output_path="\$2"
      local label="\$3"
      if [[ ! -s "\$source_path" ]]; then
        echo "ERROR: EGAPx did not produce expected \$label output: \$source_path" >&2
        echo "       Check \$PERSIST_DIR/egapx_out/nextflow logs, executor '\$EGAPX_EXECUTOR', data version '\$EGAPX_DATA_VERSION', and container '\$EGAPX_CONTAINER'." >&2
        exit 1
      fi
      cp "\$source_path" "\$output_path"
    }

    copy_required_egapx_output "\$PERSIST_DIR/egapx_out/complete.genomic.gff" egapx.complete.genomic.gff3 "complete genomic GFF"
    copy_required_egapx_output "\$PERSIST_DIR/egapx_out/complete.genomic.gtf" egapx.complete.genomic.gtf "complete genomic GTF"
    copy_required_egapx_output "\$PERSIST_DIR/egapx_out/complete.proteins.faa" egapx.complete.proteins.faa "protein FASTA"
    copy_required_egapx_output "\$PERSIST_DIR/egapx_out/complete.cds.fna" egapx.complete.cds.fna "CDS FASTA"
    copy_required_egapx_output "\$PERSIST_DIR/egapx_out/complete.transcripts.fna" egapx.complete.transcripts.fna "transcript FASTA"
    copy_required_egapx_output "\$PERSIST_DIR/egapx_out/annotated_genome.asn" egapx.annotated_genome.asn "annotated genome ASN"

    # `egapx_out` is a declared process output; when it lives in PERSIST_DIR
    # (EGAPX_WORK_DIR mode) it must also appear in the task work dir for
    # Nextflow to stage it -- a symlink avoids copying a multi-GB tree.
    if [[ "\$PERSIST_DIR" != "\$TASK_PWD" ]]; then
      ln -s "\$PERSIST_DIR/egapx_out" egapx_out
    fi

    printf '"%s":\\n  egapx: "%s"\\n  egapx_runner_revision: "%s"\\n  egapx_container: "%s"\\n  egapx_data_version: "%s"\\n' \\
      "${task.process}" "\$EGAPX_VERSION" "\$EGAPX_REVISION" "\$EGAPX_CONTAINER" "\$EGAPX_DATA_VERSION" > versions.yml
    """

  stub:
    """
    mkdir -p egapx_out/GNOMON egapx_out/stats egapx_out/validated egapx_out/nextflow
    printf "##gff-version 3\\nchr1\\tEGAPx\\tgene\\t1\\t10\\t.\\t+\\t.\\tID=egapx_stub_gene\\n" > egapx.complete.genomic.gff3
    printf "chr1\\tEGAPx\\tgene\\t1\\t10\\t.\\t+\\t.\\tgene_id \\"egapx_stub_gene\\";\\n" > egapx.complete.genomic.gtf
    printf ">egapx_stub_protein\\nM\\n" > egapx.complete.proteins.faa
    printf ">egapx_stub_cds\\nATG\\n" > egapx.complete.cds.fna
    printf ">egapx_stub_transcript\\nATG\\n" > egapx.complete.transcripts.fna
    printf "EGAPx ASN stub\\n" > egapx.annotated_genome.asn
    cp egapx.complete.genomic.gff3 egapx_out/complete.genomic.gff
    cp egapx.complete.genomic.gtf egapx_out/complete.genomic.gtf
    cp egapx.complete.proteins.faa egapx_out/complete.proteins.faa
    cp egapx.complete.cds.fna egapx_out/complete.cds.fna
    cp egapx.complete.transcripts.fna egapx_out/complete.transcripts.fna
    cp egapx.annotated_genome.asn egapx_out/annotated_genome.asn
    printf '"%s":\\n  egapx: "%s"\\n  egapx_runner_revision: "%s"\\n  egapx_container: "%s"\\n  egapx_data_version: "%s"\\n' \\
      "${task.process}" "${params.egapx_version}" "${params.egapx_revision}" "${params.container_egapx}" "${params.egapx_data_version}" > versions.yml
    """
}
