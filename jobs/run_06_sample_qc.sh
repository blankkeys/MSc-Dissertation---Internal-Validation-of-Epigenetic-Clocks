#!/bin/bash
#SBATCH --job-name=sample_qc 
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G

# Slurm job for 06/06_sample_qc.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/04_preprocessing/06_sample_qc.r


