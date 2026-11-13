#!/bin/bash

#SBATCH --partition=Centaurus
#SBATCH --job-name=a2-q1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --exclusive
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=a2-q1_%j.log

set -e

cd "$SLURM_SUBMIT_DIR"

echo "Node: $(hostname)"
echo "Date: $(date)"
echo "Compiler:"
g++ --version | head -n 1

echo "Building program..."
make clean
make

echo "Running Q1..."
cd src
./pthreads
