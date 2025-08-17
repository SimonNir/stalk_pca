#!/bin/bash

#SBATCH -N 1

#SBATCH -n 32

##SBATCH --mem-per-cpu=1G     # for most cases, default memory expectations should be fine

#SBATCH -t 8:00:00
#SBATCH -o neb.out
##SBATCH --account=brubenst-condo
#SBATCH --job-name=stalk_test
#SBATCH --mail-user=simon_nirenberg@brown.edu
#SBATCH --mail-type=END

#SBATCH --partition=batch

module load python/3.11.0s-ixrhc3q
source ~/stalk_venv/bin/activate
export PYTHONPATH=~/stalk_pca:$PYTHONPATH

python -u run1_neb.py 
