#!/bin/bash
#SBATCH --partition=physical
#SBATCH --constraint=physg5 
#SBATCH --time=2:00:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 
#SBATCH --mem=512G
#SBATCH --job-name=ldtran_test_run

module load gcc/10.1.0
g++ -fopenmp -o kseqalign_para_loops kseqalign_para_loops.cpp
./kseqalign_para_loops < mseq-big13-example.dat