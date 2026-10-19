#!/bin/bash
#SBATCH --job-name=50AM-seed-attribution
#SBATCH --account=grp_pyoung3
#SBATCH --partition=htc
#SBATCH --qos=public
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/rpraka15/50AM_full_nucleosynthesis_comparison/logs/seeds-%j.log
#SBATCH --error=/scratch/rpraka15/50AM_full_nucleosynthesis_comparison/logs/seeds-%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
cd /scratch/rpraka15/50AM_full_nucleosynthesis_comparison
/home/rpraka15/.conda/envs/ccsn/bin/python scripts/seed_attribution.py
