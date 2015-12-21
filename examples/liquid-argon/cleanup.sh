#!/bin/bash

# Script to clean up old simulation results and figures.
#
# WARNING: THIS WILL DELETE ALL SIMULATION RESULTS AND FIGURES

# Clean up compiled Python files.
rm *.pyc

# Delete data
rm -rf data.nc

# Delete figures
rm -rf figures/


