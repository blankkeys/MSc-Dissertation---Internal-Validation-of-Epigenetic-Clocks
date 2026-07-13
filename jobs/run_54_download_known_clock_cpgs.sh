#!/bin/bash
#SBATCH --job-name=download_clock_cpgs
#SBATCH --output=logs/download_clock_cpgs_%j.out
#SBATCH --error=logs/download_clock_cpgs_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Download published epigenetic clock CpG lists from pyaging
# The temporary Python environment is made inside the job scratch space

set -euo pipefail

cd /data/home/bt25127/Msc_Dissertation

mkdir -p logs
mkdir -p data/reference

venv_dir="${TMPDIR:-/tmp}/pyaging_venv_${SLURM_JOB_ID}"

python3 -m venv "$venv_dir"
source "$venv_dir/bin/activate"

python -m pip install --upgrade pip
python -m pip install --no-cache-dir pyaging

python scripts/09_exploratory_biological_interpretation/00_download_known_clock_cpgs.py
