#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N split_boot
#$ -q omni
#$ -pe sm 24
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1-10

source activate smcpp_mod

# run split
smc++ split --cores 24 -o split_${SGE_TASK_ID}/ boot_${SGE_TASK_ID}/model.final.json boot2_${SGE_TASK_ID}/model.final.json \
--timepoints 5e3 2e6 bootstrap_${SGE_TASK_ID}/*total.smc

