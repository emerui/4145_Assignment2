#!/bin/bash

#SBATCH --partition=Centaurus
#SBATCH --job-name=a2-q2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --exclusive
#SBATCH --mem=16G
#SBATCH --time=00:50:00
#SBATCH --output=a2-q2_%j.log

set -e

cd "$SLURM_SUBMIT_DIR"

echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Compiler:"
g++ --version | head -n 1

echo "Building program..."
make clean
make
echo "Running Part A (pthreads)..."
cd src
./pthreads_out

echo "Running Part B (OpenMP)..."
export OMP_NUM_THREADS=36
export OMP_PROC_BIND=true
export OMP_PLACES=cores
./openmp_out

echo "Both programs finished."
