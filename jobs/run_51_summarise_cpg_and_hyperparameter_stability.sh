#!/bin/bash
#SBATCH --job-name=cpg_stability
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Slurm job for 51_summarise_cpg_and_hyperparameter_stability.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs
mkdir -p results/analysis

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/09_validation_informed_clocks/51_summarise_cpg_and_hyperparameter_stability.r
