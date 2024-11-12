#!/bin/bash
#SBATCH --job-name=sim-study
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G		        # memory per cpu-core (4G is default)
#SBATCH --time=1-00:00:00           # Total run time limit (D-HH:MM:SS))
#SBATCH --partition=preemptable     # Partition
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jrengle@emory.edu
#SBATCH --output=HPC/Output/slurm-%A.%a.out   # Output file
#SBATCH --error=HPC/Error/slurm-%A.%a.err     # Error file
#SBATCH --array=1-200

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID."
echo "Executing on the machine:" $(hostname)

export OMP_NUM_THREADS=1

# Load R module
module load R/4.4.0

# The key is passed in as the first argument from run-all-models.sh
key=$1

# The seed is set based on the array task ID
seed=$SLURM_ARRAY_TASK_ID

# Print current key and seed being processed (for debugging/logging purposes)
echo "Processing key: $key with seed: $seed"

# Run fit-model.R for the current key and seed
Rscript R/run-sim.R $key $seed
