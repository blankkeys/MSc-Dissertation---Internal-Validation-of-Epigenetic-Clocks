#!/bin/bash
#SBATCH --job-name=make_GSE51032_sample_sheet
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

# Slurm job for 25_make_external_gse51032_sample_sheet.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/08_external_validation/08B_check_and_match_metadata/25_make_external_gse51032_sample_sheet.r
