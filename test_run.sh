#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=01-12:00:00
#SBATCH -o output_%j.out                    
#SBATCH -e error_%j.err                   
#SBATCH --job-name=test
#SBATCH --mem=200G
#SBATCH --partition=short

echo "Job started at: $(date)"
module load anaconda3/2024.06
eval "$(conda shell.bash hook)"
conda activate 2d-sweep

L1=100
L2=100
rho=100
m=0.25
s=0.005

echo "SLiM sim started at: $(date)"
slim -d L1=$L1 -d L2=$L2 -d rho=$rho -d m=$m -d s=$s 2d_spatial_sweep.slim
echo "slim sim finished at: $(date)"
python -u post_slim_process.py
echo "finished pyslim + msprime post-processing at: $(date)"
