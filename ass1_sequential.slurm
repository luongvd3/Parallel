#!/bin/bash
#SBATCH --partition=physical
#SBATCH --constraint=physg5 
#SBATCH --time=2:00:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=32G
#SBATCH --job-name=ldtran_test_run

module load gcc/10.1.0
export OMP_PROC_BIND=true
export OMP_NESTED=TRUE
g++ -fopenmp  -lnuma -Wall -O3 -o kseqalign_sequential kseqalign_sequential.cpp
./kseqalign_sequential < mseq-big13-example.dat