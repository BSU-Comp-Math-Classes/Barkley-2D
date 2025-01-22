#!/bin/bash
###
###
#SBATCH --time=00:10:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --partition=bsudfq
#SBATCH --job-name=barkley
#SBATCH --output=barkley.o%j

#module load gcc mpich slurm


./barkley 100

