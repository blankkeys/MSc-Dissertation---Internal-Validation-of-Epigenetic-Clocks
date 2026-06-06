#!/bin/bash
#SBATCH --job-name=remove_unknown_age
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

# Slurm job for 19_remove_unknown_age_samples.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/05_modelling_preparation/19_remove_unknown_age_samples.r
