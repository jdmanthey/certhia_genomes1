#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N diplo_1
#$ -q omni
#$ -pe sm 12
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -t 1-11

source activate selection

revised_number="$((${SGE_TASK_ID}-1))"

python ~/diploSHIC/diploSHIC.py fvecSim diploid albescens_hard_${revised_number}.txt.gz \
albescens_hard_${revised_number}.fvec --totalPhysLen 220000 \
--maskFileName diplo_mask.fasta.gz --chrArmsForMasking all --numSubWins 11


python ~/diploSHIC/diploSHIC.py fvecSim diploid albescens_soft_${revised_number}.txt.gz \
albescens_soft_${revised_number}.fvec --totalPhysLen 220000 \
--maskFileName diplo_mask.fasta.gz --chrArmsForMasking all --numSubWins 11


python ~/diploSHIC/diploSHIC.py fvecSim diploid americana_hard_${revised_number}.txt.gz \
americana_hard_${revised_number}.fvec --totalPhysLen 220000 \
--maskFileName diplo_mask.fasta.gz --chrArmsForMasking all --numSubWins 11


python ~/diploSHIC/diploSHIC.py fvecSim diploid americana_soft_${revised_number}.txt.gz \
americana_soft_${revised_number}.fvec --totalPhysLen 220000 \
--maskFileName diplo_mask.fasta.gz --chrArmsForMasking all --numSubWins 11

