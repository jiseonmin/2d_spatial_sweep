#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2-00:00:00
#SBATCH -o output_%j.out                    
#SBATCH -e error_%j.err                   
#SBATCH --job-name=test
#SBATCH --mem=100G
#SBATCH --partition=short

module load anaconda3/2024.06
eval "$(conda shell.bash hook)"
conda activate 2d-sweep

L1=500
L2=500
rho=500
m=0.25
s=0.05

slim -d L1=$L1 -d L2=$L2 -d rho=$rho -d m=$m -d s=$s 2d_spatial_sweep.slim
