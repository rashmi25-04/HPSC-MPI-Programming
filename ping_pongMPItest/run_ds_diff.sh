#!/bin/bash

#SBATCH --job-name ping_pong
#SBATCH --qos janus-debug

#SBATCH --nodes 2
#SBATCH --ntasks-per-node 1
#SBATCH --nodelist=node[0419],node[1630]
#SBATCH --time 00:20:00
#SBATCH --output=hw3-diffswitch-%j.out
#SBATCH --error=hw3-diffswitch-%j.err

# the slurm module provides the srun command
module load slurm

module load openmpi


srun ping_pong
