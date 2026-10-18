#!/bin/bash
#SBATCH --job-name=compare-isotope-tables
#SBATCH --account=grp_pyoung3
#SBATCH --partition=htc
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/rpraka15/file_comparison_20261018/slurm-%j.log
#SBATCH --error=/scratch/rpraka15/file_comparison_20261018/slurm-%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
exec /home/rpraka15/.conda/envs/ccsn/bin/python -u /scratch/rpraka15/file_comparison_20261018/compare_outputs.py \
  --final /scratch/rpraka15/final_chunks.out \
  --plotting /scratch/rpraka15/50Am_plotting_all_isotopes.out \
  --output "/scratch/rpraka15/file_comparison_20261018/results-${SLURM_JOB_ID}" \
  --expected-final-bytes 1438981266 \
  --expected-final-sha256 3c9409582bd909f119c3484b9ae7ea42895aed6da5993df0ffdd974775fb78ca \
  --expected-plotting-sha256 795949280886f028bd9391f8f43af6a1d38e278014104e11168cd4b10ddc0399
