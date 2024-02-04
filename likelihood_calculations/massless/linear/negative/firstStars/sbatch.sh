#!/bin/bash
# 
# CompecTA (c) 2017
#
# You should only work under the /scratch/users/<username> directory.
#
# Example job submission script
#
# TODO:
#   - Set name of the job below changing "Test" value.
#   - Set the requested number of tasks (cpu cores) with --ntasks parameter.
#   - Select the partition (queue) you want to run the job in:
#     - short : For jobs that have maximum run time of 120 mins. Has higher priority.
#     - long  : For jobs that have maximum run time of 7 days. Lower priority than short.
#     - longer: For testing purposes, queue has 31 days limit but only 3 nodes.
#   - Set the required time limit for the job with --time parameter.
#     - Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds"
#   - Put this script and all the input file under the same directory.
#   - Set the required parameters, input and output file names below.
#   - If you do not want mail please remove the line that has --mail-type
#   - Put this script and all the input file under the same directory.
#   - Submit this file using:
#      sbatch examle_submit.sh

# -= Resources =-
#
#SBATCH --job-name=firststar3
#SBATCH --partition=long
#SBATCH --output=./first3
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --time=4-00:00:00

################################################################################
##################### !!! DO NOT EDIT BELOW THIS LINE !!! ######################
################################################################################
ulimit -s unlimited
ulimit -l unlimited
ulimit -a
# Put compiled binary command below
cat /truba/home/edemirboga/vectortensor/apc524_idsolve/IDsolve.fparam ./rtparam > ./param
cp /truba/home/edemirboga/vectortensor/apc524_idsolve/FindScalarized FindScalarized
module load centos7.9/comp/gcc/7
./FindScalarized param


