#!/bin/bash
#SBATCH --job-name=fit_validation_clocks
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G

# Slurm job for 33_fit_validation_informed_clocks.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs
mkdir -p results/modelling/validation_informed_clocks

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/09_analysis/33_fit_validation_informed_clocks.r
