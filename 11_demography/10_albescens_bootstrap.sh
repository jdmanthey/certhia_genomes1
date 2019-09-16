#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N alb_boot
#$ -q omni
#$ -pe sm 24
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1-10

source activate smcpp_mod

# run with albescens
smc++ cv 2.506e-09 -o boot_${SGE_TASK_ID} --cores 24 --knots 16 --timepoints 5e3 2e6 --regularization-penalty 6 \
--spline cubic bootstrap_${SGE_TASK_ID}/*albescens.smc
