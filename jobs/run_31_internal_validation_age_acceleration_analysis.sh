#!/bin/bash
#SBATCH --job-name=age_acceleration_analysis
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Slurm job for 31_internal_validation_age_acceleration_analysis.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs
mkdir -p results/analysis

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/09_analysis/31_internal_validation_age_acceleration_analysis.r
