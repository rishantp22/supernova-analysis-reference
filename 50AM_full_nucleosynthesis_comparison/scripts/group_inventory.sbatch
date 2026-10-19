#!/bin/bash
#SBATCH --job-name=50AM-all-group-files
#SBATCH --account=grp_pyoung3
#SBATCH --partition=htc
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/rpraka15/50AM_full_nucleosynthesis_comparison/logs/group-inventory-%j.log
#SBATCH --error=/scratch/rpraka15/50AM_full_nucleosynthesis_comparison/logs/group-inventory-%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /scratch/rpraka15/50AM_full_nucleosynthesis_comparison
exec /home/rpraka15/.conda/envs/ccsn/bin/python -u scripts/scan_all_accessible.py
