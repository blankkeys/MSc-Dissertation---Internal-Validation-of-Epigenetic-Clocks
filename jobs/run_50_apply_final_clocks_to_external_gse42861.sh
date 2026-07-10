#!/bin/bash
#SBATCH --job-name=compare_validation_clocks
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G

# Slurm job for 50_apply_final_clocks_to_external_gse42861.r

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs
mkdir -p results/external_validation/validation_informed_clocks

apptainer exec --cleanenv containers/bioconductor_3_22.sif \
  Rscript scripts/09_validation_informed_clocks/50_apply_final_clocks_to_external_gse42861.r
