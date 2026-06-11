#!/bin/bash
#SBATCH --job-name=remove_failed_GSE42861
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G

# Slurm job for 36_external_gse42861_remove_failed_samples.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/08_external_validation/08D_preprocessing/36_external_gse42861_remove_failed_samples.r
