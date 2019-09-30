#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N fvec_ame1
#$ -q omni
#$ -pe sm 4
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1-11

source activate selection2

modnum=$((SGE_TASK_ID-1))

python ~/diploSHIC/diploSHIC.py fvecSim diploid americana_soft_${modnum}.txt.gz \
americana_soft_${modnum}.fvec --totalPhysLen 220000 \
--maskFileName diplo_mask.fasta.gz --chrArmsForMasking all --numSubWins 11

python ~/diploSHIC/diploSHIC.py fvecSim diploid americana_hard_${modnum}.txt.gz \
americana_hard_${modnum}.fvec --totalPhysLen 220000 \
--maskFileName diplo_mask.fasta.gz --chrArmsForMasking all --numSubWins 11

