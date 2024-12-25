#!/bin/bash

 

#SBATCH --partition=highmem
##SBATCH --exclusive
#SBATCH -n 64
#SBATCH --ntasks-per-node=64

#SBATCH --mem=1000000

 

source ~/.bashrc
module load icc/20.4
module load openmpi/4.1.2

 

export OMPI_NUM_THREADS=1
mpirun ./a.out
