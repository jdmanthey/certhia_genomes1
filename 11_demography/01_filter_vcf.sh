#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N filter
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=12G
#$ -t 1-25

vcftools --vcf ${SGE_TASK_ID}.rn.vcf --max-missing 0.5 --minQ 20 --minGQ 20 --minDP 5 \
--max-meanDP 60 --recode --remove-indels --keep smc_keep.txt --out ${SGE_TASK_ID}
