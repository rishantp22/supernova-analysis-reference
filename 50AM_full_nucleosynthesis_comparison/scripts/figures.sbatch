#!/bin/bash
#SBATCH --job-name=50AM-audit-figures
#SBATCH --account=grp_pyoung3
#SBATCH --partition=htc
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/rpraka15/50AM_full_nucleosynthesis_comparison/logs/figures-%j.log
#SBATCH --error=/scratch/rpraka15/50AM_full_nucleosynthesis_comparison/logs/figures-%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=/scratch/rpraka15/50AM_full_nucleosynthesis_comparison/figures/.mplconfig
cd /scratch/rpraka15/50AM_full_nucleosynthesis_comparison
exec /home/rpraka15/.conda/envs/ccsn/bin/python -u scripts/figures.py
