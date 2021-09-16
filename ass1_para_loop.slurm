#!/bin/bash
#SBATCH --partition=physical
#SBATCH --constraint=physg5 
#SBATCH --time=0:30:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --job-name=ldtran_test_run

module load gcc/10.1.0
export OMP_PROC_BIND=true
export OMP_NESTED=TRUE
export OMP_PLACES={0:2},{2:2},{4:2},{6:2}
g++ -fopenmp  -lnuma -Wall -O3 -o kseqalign_para_loops kseqalign_para_loops.cpp
./kseqalign_para_loops < mseq-big13-example.dat