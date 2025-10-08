#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00-12:00:00
#SBATCH -o output_%j.out                    
#SBATCH -e error_%j.err                   
#SBATCH --job-name=test
#SBATCH --mem=50G
#SBATCH --partition=short

echo "Job started at: $(date)"
module load anaconda3/2024.06
eval "$(conda shell.bash hook)"
conda activate 2d-sweep

python -u test_big_merge.py
echo "finished testing merger: $(date)"
