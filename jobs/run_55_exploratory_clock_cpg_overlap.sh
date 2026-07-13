#!/bin/bash
#SBATCH --job-name=cpg_overlap
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

# Slurm job for exploratory CpG overlap with established epigenetic clocks

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs
mkdir -p data/reference
mkdir -p results/exploratory_cpg_overlap

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/09_exploratory_biological_interpretation/01_clock_cpg_overlap.R
