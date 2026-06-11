#!/bin/bash
#SBATCH --job-name=calibration_bias
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Slurm job for 29_internal_validation_calibration_and_bias.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/07_internal_validation/29_internal_validation_calibration_and_bias.r
