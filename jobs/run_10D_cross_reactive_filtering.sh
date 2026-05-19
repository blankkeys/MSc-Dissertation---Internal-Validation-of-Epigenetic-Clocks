#!/bin/bash
#SBATCH --job-name=cross_reactive_filtering
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G

# Slurm job for 10D_cross_reactive_filtering.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/04_preprocessing/10D_cross_reactive_filtering.r
