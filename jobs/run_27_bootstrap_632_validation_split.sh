#!/bin/bash
#SBATCH --job-name=bootstrap_632_split
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --array=1-10

# Slurm array job for 27_bootstrap_632_validation.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs
mkdir -p results/internal_validation

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/07_internal_validation/27_bootstrap_632_validation.r
