#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N ca_raxml
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe sm 1
#$ -P quanah
#$ -l h_rt=01:00:00
#$ -l h_vmem=8G
#$ -t 1-18470

raxmlHPC-PTHREADS-SSE3 -T 1 -f a -x 50 -m GTRGAMMA -p 253 -N 100 \
-s /lustre/scratch/jmanthey/01_certhia_genomics/certhia_fasta/${SGE_TASK_ID}.fasta \
-n certhia_${SGE_TASK_ID}.tre -w /lustre/scratch/jmanthey/01_certhia_genomics/certhia_fasta/

rm RAxML_bestTree.certhia_${SGE_TASK_ID}.tre

rm RAxML_bipartitionsBranchLabels.certhia_${SGE_TASK_ID}.tre

rm RAxML_bootstrap.certhia_${SGE_TASK_ID}.tre

rm RAxML_info.certhia_${SGE_TASK_ID}.tre
