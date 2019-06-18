#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N smc_ame2
#$ -q omni
#$ -pe sm 24
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G

source activate smcpp_mod

# run with americana
smc++ cv 2.506e-09 -o americana_2 --cores 24 --knots 16 --timepoints 5e3 2e6 --regularization-penalty 6 \
--spline cubic *americana.smc
