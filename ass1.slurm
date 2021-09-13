#!/bin/bash
#SBATCH --partition=physical
#SBATCH --constraint=physg5 
#SBATCH --time=0:05:00 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem=32G
#SBATCH --job-name= ldtran_test_run

module load gcc/10.1.0
./kseqalign < mseq-big13-example.dat