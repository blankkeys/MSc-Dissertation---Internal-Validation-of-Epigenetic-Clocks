#!/bin/bash
#SBATCH --job-name=nested_validation_split
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --array=1-10

# Slurm array job for 26_nested_cross_validation.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs
mkdir -p results/internal_validation

apptainer exec --cleanenv \
  --env SLURM_ARRAY_TASK_ID="${SLURM_ARRAY_TASK_ID}" \
  containers/bioconductor_3_22.sif \
  Rscript scripts/07_internal_validation/26_nested_cross_validation.r
