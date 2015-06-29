#!/bin/bash

# Script to clean up old simulation results and figures.
#
# WARNING: THIS WILL DELETE ALL SIMULATION RESULTS AND FIGURES

# Delete data
rm -rf data/

# Delete figures
rm -rf figures/

# Remove conda environment
rm -rf conda-env/


