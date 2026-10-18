#!/bin/bash
#SBATCH --job-name=compare-isotope-tables
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.log
#SBATCH --error=slurm-%j.err
set -euo pipefail
if [ "$#" -ne 3 ]; then
  echo 'Usage: sbatch slurm/compare_tables.sbatch FINAL_TABLE PLOTTING_TABLE NEW_OUTPUT_DIR' >&2
  exit 2
fi
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
exec "${PYTHON_BIN:-python3}" -u comparison/compare_outputs.py \
  --final "$1" --plotting "$2" --output "$3"
