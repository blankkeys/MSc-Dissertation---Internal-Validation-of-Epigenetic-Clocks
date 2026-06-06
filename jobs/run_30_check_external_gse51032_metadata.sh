#!/bin/bash
#SBATCH --job-name=check_GSE51032_metadata
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Slurm job for 30_check_external_gse51032_metadata.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs
mkdir -p results/external_validation

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/08_external_validation/08B_check_and_match_metadata/30_check_external_gse51032_metadata.r
