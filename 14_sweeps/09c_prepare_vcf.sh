#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N diplo_vcf
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1-26

vcftools --vcf ${SGE_TASK_ID}_final.recode.vcf --max-missing 0.5 --keep americana.txt \
--min-alleles 2 --max-alleles 2 --mac 1 --max-meanDP 60 --recode --recode-INFO-all --out ${SGE_TASK_ID}_americana

vcftools --vcf ${SGE_TASK_ID}_diplo.recode.vcf --max-missing 0.5 --keep albescens.txt \
--min-alleles 2 --max-alleles 2 --mac 1 --max-meanDP 60 --recode --recode-INFO-all --out ${SGE_TASK_ID}_albescens
