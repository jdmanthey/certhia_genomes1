#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N prep_smc
#$ -q omni
#$ -pe sm 2
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=12G
#$ -t 1-25

source activate smcpp_mod

j=$(gzip -cd ${SGE_TASK_ID}_mask.bed.gz | head -n1 | cut -f1)

smc++ vcf2smc -d 22 22 ${SGE_TASK_ID}.recode.vcf.gz ${SGE_TASK_ID}_a_americana.smc $j \
americana:13,14,15,16,17,18,19,20,21,22,23,24 -m ${SGE_TASK_ID}_mask.bed.gz;

smc++ vcf2smc -d 23 23 ${SGE_TASK_ID}.recode.vcf.gz ${SGE_TASK_ID}_b_americana.smc $j \
americana:13,14,15,16,17,18,19,20,21,22,23,24 -m ${SGE_TASK_ID}_mask.bed.gz;

smc++ vcf2smc -d 10 10 ${SGE_TASK_ID}.recode.vcf.gz ${SGE_TASK_ID}_a_albescens.smc $j \
albescens:2,7,8,9,10,11,12 -m ${SGE_TASK_ID}_mask.bed.gz;

smc++ vcf2smc -d 11 11 ${SGE_TASK_ID}.recode.vcf.gz ${SGE_TASK_ID}_b_albescens.smc $j \
albescens:2,7,8,9,10,11,12 -m ${SGE_TASK_ID}_mask.bed.gz;

smc++ vcf2smc ${SGE_TASK_ID}.recode.vcf.gz ${SGE_TASK_ID}_a_total.smc $j \
albescens:2,7,8,9,10,11,12 americana:13,14,15,16,17,18,19,20,21,22,23,24 -m ${SGE_TASK_ID}_mask.bed.gz;

smc++ vcf2smc ${SGE_TASK_ID}.recode.vcf.gz ${SGE_TASK_ID}_b_total.smc $j \
americana:13,14,15,16,17,18,19,20,21,22,23,24 albescens:2,7,8,9,10,11,12 -m ${SGE_TASK_ID}_mask.bed.gz;

