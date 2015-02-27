#!/bin/bash
#  Batch script for submitting ocores to MSKCC Torque/Moab cluster.
#
#
# walltime : maximum wall clock time (hh:mm:ss)
#PBS -l walltime=12:00:00
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
#
# filename for standard output (default = <job_name>.o<job_id>)
# at end of job, it is in directory from which qsub was executed
# remove extra ## from the line below if you want to name your own file
##PBS -o /cbio/jclab/home/chodera/ocores/output

cd "$PBS_O_WORKDIR"

export CUDA_VISIBLE_DEVICES=`cat $PBS_GPUFILE | awk -F"-gpu" '{ printf A$2;A=","}'`

date
hostname
echo CUDA_VISIBLE_DEVICES $CUDA_VISIBLE_DEVICES

python simulate.py

date

