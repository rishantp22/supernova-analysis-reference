#!/bin/bash
#SBATCH --job-name=50AM-final-validation
#SBATCH --account=grp_pyoung3
#SBATCH --partition=htc
#SBATCH --qos=public
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/rpraka15/50AM_full_nucleosynthesis_comparison/logs/final-validation-%j.log
#SBATCH --error=/scratch/rpraka15/50AM_full_nucleosynthesis_comparison/logs/final-validation-%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
cd /scratch/rpraka15/50AM_full_nucleosynthesis_comparison
/home/rpraka15/.conda/envs/ccsn/bin/python scripts/final_validation.py
