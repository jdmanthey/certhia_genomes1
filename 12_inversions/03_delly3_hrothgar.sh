#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N C_delly3
#$ -q Chewie
#$ -pe sm 4
#$ -P communitycluster
#$ -l h_rt=48:00:00
#$ -l h_vmem=10G
#$ -t 1:25


/home/jmanthey/delly_v0.8.1_linux_x86_64bit call -g /home/jmanthey/references/06_certhia_reordered.fasta \
-v /lustre/scratch/jmanthey/01_certhia_genomics/09_delly/certhia_sites.bcf \
-o /lustre/scratch/jmanthey/01_certhia_genomics/09_delly/${SGE_TASK_ID}_SV.geno.bcf \
-t INV -q 20 -s 15 \
-x /lustre/scratch/jmanthey/01_certhia_genomics/09_delly/exclusion_list.bed \
/lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}_final.bam
