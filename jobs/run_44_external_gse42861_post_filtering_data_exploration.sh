#!/bin/bash
#SBATCH --job-name=post_filter_GSE42861
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G

# Slurm job for 44_external_gse42861_post_filtering_data_exploration.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs
mkdir -p results/external_validation/qc

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/08_external_validation/08E_probe_filtering/44_external_gse42861_post_filtering_data_exploration.r
