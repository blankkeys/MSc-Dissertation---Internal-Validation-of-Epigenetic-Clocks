#!/bin/bash
#SBATCH --job-name=cross_reactive_GSE51032
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G

# Slurm job for 40_external_gse51032_cross_reactive_filtering.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs
mkdir -p results/external_validation/qc

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/08_external_validation/08E_probe_filtering/40_external_gse51032_cross_reactive_filtering.r
