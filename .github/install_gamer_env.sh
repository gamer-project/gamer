#!/bin/bash

# Install GAMER Environment for GitHub-hosted GitHub Actions runners
# This script installs dependencies needed to compile and run GAMER

set -e

echo "Installing GAMER dependencies..."

# Install OpenMPI, HDF5, and CUDA toolkit
sudo apt-get install -y openmpi-bin libopenmpi-dev libhdf5-dev nvidia-cuda-toolkit

echo ""
echo "Verifying installations..."
echo "MPI:"
which mpirun
which mpicxx
mpicxx --version

echo ""
echo "CUDA:"
which nvcc
nvcc --version

echo ""
echo "HDF5:"
h5cc -show

echo ""
echo "Configuring GAMER machine settings..."
cd "$(git rev-parse --show-toplevel)/tool/config"
bash set_settings.sh --local --machine=github_action

echo ""
echo "Verifying configuration..."
cat "$(git rev-parse --show-toplevel)/src/.local_settings"

echo ""
echo "GAMER environment setup complete!"
