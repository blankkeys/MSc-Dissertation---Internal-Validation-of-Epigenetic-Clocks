#!/bin/bash
#SBATCH --job-name=extract_beta_matrix
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G

# Slurm job for 11/11_extract_beta_matrix.r

set -euo pipefail # Exit on error, treat unset variables as errors, and fail if any command in a pipeline fails.

# Run from the project root so paths like data/GSE87571 work correctly.
cd /data/home/bt25127/Msc_Dissertation

# Make sure Slurm has somewhere to write the output and error logs.
mkdir -p logs

# Run the R script inside the Bioconductor container to ensure package compatibility.
apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/04_preprocessing/11_extract_beta_matrix.r

