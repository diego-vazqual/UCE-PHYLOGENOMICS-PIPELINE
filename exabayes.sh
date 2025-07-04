#!/bin/bash
#SBATCH --job-name=BI_lowers # Name of job
#SBATCH --output=BI_lowers_%j.out # Output
#SBATCH --error=BI_lowers_%j.err # Possible errors
#SBATCH --ntasks=16 # Number of parallel task you want to execute
#SBATCH --cpus-per-task=5 # Number of threads per task
#SBATCH --nodes=1 # How many nodes you want to use
#SBATCH --nodelist=node3 # Node you want to use
#SBATCH --mem=800G
#SBATCH --time=40-00:00:00 # Time limit for execution
#SBATCH --partition=irbio01

# Modules to load
. /etc/profile
module load anaconda/2024.02
module load exabayes/1.5.1_ompi

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Variables
ALIGNMENT="mafft-clean-nexus-internal-trimmed-gblocks-clean-50p-IQTree.phylip"
CONFIG="config.nex"

# ExaBayes command
for i in {1..4}; do
  echo "=== Executing run$i ==="
  SEED=$((1234 + i))
  mpirun -np 4 exabayes -f "$ALIGNMENT" -m DNA -c "$CONFIG" \
         -n run$i -s "$SEED" -R 1 -C 4 -T 5 -M 1 &
done

wait
