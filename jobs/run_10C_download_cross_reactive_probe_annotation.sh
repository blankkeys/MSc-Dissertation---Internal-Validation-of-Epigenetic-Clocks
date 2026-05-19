#!/bin/bash
#SBATCH --job-name=download_probe_annotation
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Slurm job for 10C_download_cross_reactive_probe_annotation.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/04_preprocessing/10C_download_cross_reactive_probe_annotation.r
