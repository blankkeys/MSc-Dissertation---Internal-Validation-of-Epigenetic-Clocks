#!/bin/bash
#SBATCH --job-name=elastic_net_regression
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G

# Slurm job for 14_elastic_net_regression.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/06_elastic_net_modelling/14_elastic_net_regression.r
