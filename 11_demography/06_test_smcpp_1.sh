#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N test_smc
#$ -q omni
#$ -pe sm 12
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

source activate smcpp_mod

# run with americana
smc++ cv 2.506e-09 -o americana_test1 --cores 12 --knots 16 --timepoints 5e2 3e6 --regularization-penalty 10 \
--spline cubic *americana.smc
