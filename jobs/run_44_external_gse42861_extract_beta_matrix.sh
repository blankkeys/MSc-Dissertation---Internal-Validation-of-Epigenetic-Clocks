#!/bin/bash
#SBATCH --job-name=beta_matrix_GSE42861
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G

# Slurm job for 44_external_gse42861_extract_beta_matrix.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs
mkdir -p results/external_validation/qc

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/08_external_validation/08F_extract_beta_matrix/44_external_gse42861_extract_beta_matrix.r
