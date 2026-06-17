#!/bin/bash
#SBATCH --job-name=validation_plots
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Slurm job for internal_validation_prediction_plots.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/descriptive_plots/internal_validation_prediction_plots.r
