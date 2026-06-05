#!/bin/bash
#SBATCH --job-name=download_GSE51032
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

# Slurm job for 22_download_external_gse51032.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/08_external_validation/08A_download_gse51032/22_download_external_gse51032.r
