#!/bin/bash

#SBATCH --job-name fox
#SBATCH --qos janus-debug

#SBATCH --nodes 64

#SBATCH --time 00:20:00
#SBATCH --output=hw3-fox-%j.out
#SBATCH --error=hw3-fox-%j.err

# the slurm module provides the srun command
module load slurm

module load openmpi


srun fox
