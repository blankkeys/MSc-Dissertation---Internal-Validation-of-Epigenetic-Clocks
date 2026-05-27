#!/bin/bash
#SBATCH --job-name=repeated_train_test
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G

# Slurm job for 17_repeated_train_test_split.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs
mkdir -p results/internal_validation

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/07_internal_validation/17_repeated_train_test_split.r
