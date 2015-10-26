#!/bin/bash

#SBATCH --job-name mm_mpi
#SBATCH --qos janus-debug

#SBATCH --nodes 144

#SBATCH --time 00:20:00
#SBATCH --output=hw3-allgather-%j.out
#SBATCH --error=hw3-allgather-%j.err

# the slurm module provides the srun command
module load slurm

module load openmpi


srun mm_mpi
