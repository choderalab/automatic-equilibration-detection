#!/bin/bash

# Script to reproduce all results of paper from scratch.
# Note that you must have called `cleanup.sh` script first to clean up old files.
#
# PREREQUISITES
#
# This requires the conda Pythjon package environment management tool.

# Create conda environment with necessary tools.
if [ ! -d conda-env ]; then
    echo "Creating new conda environment ./conda-env containing versions of tools needed for reproducibility..."
    conda config --add channels http://conda.binstar.org/omnia
    conda create --yes --quiet -p conda-env python=2.7 openmmtools=0.7.0 openmm=6.2 matplotlib=1.4 pymbar=3.0.0.beta2 netCDF4 seaborn matplotlib
fi
source activate ./conda-env

# Run simulations.
python simulate.py

# Analyze simulation data to generate figures.
if [ ! -e figures ]; then
    mkdir figures
fi
python analyze-1.py
python analyze-2.py

# Deactivate conda environment.
source deactivate


