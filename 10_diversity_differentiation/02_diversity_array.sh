#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N certhia_div
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=01:00:00
#$ -l h_vmem=8G
#$ -t 2-9836

module load intel R
Rscript div_diff_${SGE_TASK_ID}.r
