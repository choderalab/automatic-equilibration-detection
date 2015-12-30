#!/bin/bash

# Script to reproduce all results of paper from scratch.
# Note that you must have called `cleanup.sh` script first to clean up old files.
#
# PREREQUISITES
#
# This requires the conda Python package environment management tool.

# Create conda environment with necessary tools.
CONDA_ENV="argon"
if [ "$(conda env list | grep $CONDA_ENV | wc -l)" -eq "0" ]; then
    echo "Creating new conda environment '$CONDA_ENV' containing versions of tools needed for reproducibility..."
    conda config --add channels http://conda.binstar.org/omnia
    conda env create --quiet -n $CONDA_ENV --file environment.yml
fi
source activate $CONDA_ENV

# Run simulations.
echo "Running simulations..."
python simulate.py

# Analyze simulation data to generate figures.
if [ ! -e figures ]; then
    mkdir figures
fi
echo "Generating figures..."
python analyze-1.py
python analyze-2.py
python analyze-3.py
python analyze-4.py

# Deactivate conda environment.
source deactivate

echo "Finished."
echo ""
echo "You can now safely remove the conda environment with:"
echo ""
echo "conda env remove -n $CONDA_ENV"
echo ""
echo "If you have no further need of it."



