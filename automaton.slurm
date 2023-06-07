#!/bin/bash

# This script is configured for rapid turnaround in the short queue
# This only works for jobs *less than 20 minutes on 2 nodes or fewer*.
# For more nodes set "--qos=standard" and delete the reservation line.

#SBATCH --account=<account>
#SBATCH --job-name=automaton
#SBATCH --time=00:01:00
#SBATCH --nodes=2
#SBATCH --ntasks=64
#SBATCH --qos=short
#SBATCH --reservation=shortqos

# You can probably leave these options mostly as they are

#SBATCH --output=%x-%j.out
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --partition=standard

module load intel-compilers-19
module load mpt

# Launch the parallel job

srun --unbuffered --distribution=block:block --hint=nomultithread ./automaton 1234
