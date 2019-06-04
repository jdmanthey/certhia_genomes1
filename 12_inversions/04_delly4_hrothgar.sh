#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N C_delly4
#$ -q Chewie
#$ -pe sm 4
#$ -P communitycluster
#$ -l h_rt=48:00:00
#$ -l h_vmem=10G
#$ -t 1:1

bcftools merge -m id -O v \
-o /lustre/scratch/jmanthey/01_certhia_genomics/09_delly/certhia_inversions.vcf \
/lustre/scratch/jmanthey/01_certhia_genomics/09_delly/*_SV.geno.bcf

