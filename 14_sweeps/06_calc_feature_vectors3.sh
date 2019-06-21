#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N fvec_ame2
#$ -q omni
#$ -pe sm 4
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1-1

source activate selection2

python ~/diploSHIC/diploSHIC.py fvecSim diploid americana_neutral.txt.gz \
americana_neutral.fvec --totalPhysLen 110000 \
--maskFileName diplo_mask.fasta.gz --chrArmsForMasking all --numSubWins 11
