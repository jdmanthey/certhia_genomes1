#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N C_delly1
#$ -q Chewie
#$ -pe sm 4
#$ -P communitycluster
#$ -l h_rt=48:00:00
#$ -l h_vmem=10G
#$ -t 1:1

/home/jmanthey/delly_v0.8.1_linux_x86_64bit merge -n 100000000 -m 1000 -n 1000000 \
-o /lustre/scratch/jmanthey/01_certhia_genomics/09_delly/certhia_sites.bcf \
/lustre/scratch/jmanthey/01_certhia_genomics/09_delly/*_SV.bcf
