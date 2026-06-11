#!/bin/bash
#SBATCH --job-name=import_GSE42861
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G

# Slurm job for 34_import_external_gse42861_with_minfi.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/08_external_validation/08C_import_data_into_R_with_minfi/34_import_external_gse42861_with_minfi.r
