#!/bin/bash

#SBATCH --job-name ping_pong
#SBATCH --qos janus-debug

#SBATCH --nodes 2
#SBATCH --ntasks-per-node 1
#SBATCH --nodelist=node[0207],node[0209]
#SBATCH --time 00:20:00
#SBATCH --output=hw3-sameswitch-%j.out
#SBATCH --error=hw3-sameswitch-%j.err

# the slurm module provides the srun command
module load slurm

module load openmpi


srun ping_pong
