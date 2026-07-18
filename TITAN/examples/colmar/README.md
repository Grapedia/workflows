# Colmar Slurm/Apptainer Example

This directory keeps the site-specific Colmar production launcher and Slurm
configuration outside the generic TITAN entry points.

Files:

* `launch_TITAN_serveur_colmar.sh`: sbatch launcher used on the Colmar cluster.
* `slurm_apptainer.config`: Colmar production profile overlay with local node,
  cache and resource settings.

These files are examples, not portable defaults. For a new site, copy them to a
local path, adjust the Slurm partition/node policy, input paths, cache paths and
enabled optional branches, then launch TITAN from the repository root.

The active local copies used by production runs remain ignored at:

* `launch_TITAN_serveur_colmar.sh`
* `data/slurm_apptainer.config`
