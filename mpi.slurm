#!/bin/bash
#SBATCH --error=err-%j.err
#SBATCH --output=snowy-8-16-%j.out
#SBATCH --partition=snowy
#SBATCH --time=0:4:00
#SBATCH --nodes=8
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32G
#SBATCH --job-name=pure_mpi

# You need to load a compiler before openmpi.
module load gcc/8.3.0

module load openmpi/3.1.4

echo "Scheduler 8 nodes 1 task 16 cpus (snowy)"
mpicxx -o binary openmpi.cpp -fopenmp -O3
mpirun -n 8 binary < data1.dat


