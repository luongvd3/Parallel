```sh
#!/bin/bash
#SBATCH --partition=physical<- Where do you want it to run? [snowy,physical]
#SBATCH --constraint=physg5 <- What specific machine?#SBATCH --time=0:05:00  <----- How long should it run before we kill it
#SBATCH --ntasks=1 <---------- How many programs to run at once? (1 for openmp)
#SBATCH --cpus-per-task=24 <-- How many CPUs per program?
#SBATCH --mem=32G <----------- How much memory per task?
#SBATCH --job-name=jobname <-- What name do you wanna give it?

module load gcc/10.1.0 # <---- Load the compiler module we need
g++ -fopen -o kseqalign kseqalign.cpp < mseq-big13-example.dat #<- Compile the file
./kseqalign # <---------------- Exec the file