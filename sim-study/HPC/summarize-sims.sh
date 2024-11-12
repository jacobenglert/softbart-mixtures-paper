#!/bin/bash
#SBATCH --job-name=sum-sim                # create a name for the job
#SBATCH --nodes=1                         # Node count (number of computers)
#SBATCH --ntasks-per-node=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1                 # cpu-cores per task (>1 if multi-threaded tasks, internal parallelization within R for example)
#SBATCH --mem-per-cpu=10G                  # memory per cpu-core (4G is default)
#SBATCH --mail-type=ALL                   # Email received when job fails only (can also set as start, etc.)
#SBATCH --mail-user=jrengle@emory.edu     # Email to receive
#SBATCH --output=HPC/Output/sum-sim.out   # Output file
#SBATCH --error=HPC/Error/sum-sim.err     # Error file

export OMP_NUM_THREADS=1

module load R/4.4.0

Rscript R/summarize-sims.R
