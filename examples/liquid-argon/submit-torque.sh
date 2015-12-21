#!/bin/bash
#  Batch script for Torque/Moab cluster.
#
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=24:00:00
#
# join stdout and stderr
#PBS -j oe
#
# spool output immediately
##PBS -k oe
#
# specify queue
#PBS -q gpu
#
# nodes: number of 32-hyperthread nodes
#   ppn: how many cores per node to use (1 through 32)
#       (you are always charged for the entire node)
#PBS -l nodes=1:ppn=1:gpus=1:shared
#
# export all my environment variables to the job
##PBS -V
#
# job name (default = name of script file)
#PBS -N argon

cd "$PBS_O_WORKDIR"

date
hostname

./cleanup.sh
./reproduce.sh

date

