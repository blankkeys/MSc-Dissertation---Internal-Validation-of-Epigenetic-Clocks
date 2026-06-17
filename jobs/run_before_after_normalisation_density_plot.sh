#!/bin/bash
#SBATCH --job-name=normalisation_density
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Slurm job for before_after_normalisation_density_plot.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/descriptive_plots/before_after_normalisation_density_plot.r
