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

gzip ${SGE_TASK_ID}_albescens.recode.vcf

gzip ${SGE_TASK_ID}_americana.recode.vcf
