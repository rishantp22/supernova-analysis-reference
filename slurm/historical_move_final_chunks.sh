#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --partition=htc
#SBATCH --qos=public
#SBATCH --mem=50G
#SBATCH -o slurm_log/slurm.%A_%a.out
#SBATCH -o slurm_log/slurm.%A_%a.err

# Change to the source directory
cd /scratch/svierlin/prism-v02 || exit 1

# Copy final_chunks recursively to destination
cp -r final_chunks /scratch/rpraka15/

# Print completion message
echo "done"
