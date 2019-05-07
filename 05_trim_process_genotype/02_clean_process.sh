#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N certhia1
#$ -q omni
#$ -pe sm 4
#$ -P quanah
#$ -l h_rt=48:00:00
#$ -l h_vmem=12G
#$ -t 1-25

module load intel java bwa samtools

/lustre/work/jmanthey/bbmap/bbduk.sh in1=/lustre/scratch/jmanthey/01_certhia_genomics/00_fastq/${SGE_TASK_ID}_R1.fastq.gz in2=/lustre/scratch/jmanthey/01_certhia_genomics/00_fastq/${SGE_TASK_ID}_R2.fastq.gz out1=/lustre/scratch/jmanthey/01_certhia_genomics/01_cleaned/${SGE_TASK_ID}_R1.fastq.gz out2=/lustre/scratch/jmanthey/01_certhia_genomics/01_cleaned/${SGE_TASK_ID}_R2.fastq.gz minlen=50 ftl=10 qtrim=rl trimq=10 ktrim=r k=25 mink=7 ref=/lustre/work/jmanthey/bbmap/resources/adapters.fa hdist=1 tbo tpe

bwa mem -t 4 /home/jmanthey/references/06_certhia_reordered.fasta /lustre/scratch/jmanthey/01_certhia_genomics/01_cleaned/${SGE_TASK_ID}_R1.fastq.gz /lustre/scratch/jmanthey/01_certhia_genomics/01_cleaned/${SGE_TASK_ID}_R2.fastq.gz > /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}.sam

samtools view -b -S -o /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}.bam /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}.sam

rm /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}.sam

/lustre/work/jmanthey/gatk-4.1.0.0/gatk CleanSam -I /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}.bam -O /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}_cleaned.bam

rm /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}.bam

/lustre/work/jmanthey/gatk-4.1.0.0/gatk SortSam -I /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}_cleaned.bam -O /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}_cleaned_sorted.bam --SORT_ORDER coordinate

rm /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}_cleaned.bam

/lustre/work/jmanthey/gatk-4.1.0.0/gatk AddOrReplaceReadGroups -I /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}_cleaned_sorted.bam -O /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}_cleaned_sorted_rg.bam --RGLB 1 --RGPL illumina --RGPU unit1 --RGSM ${SGE_TASK_ID}

rm /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}_cleaned_sorted.bam

/lustre/work/jmanthey/gatk-4.1.0.0/gatk MarkDuplicates --REMOVE_DUPLICATES true --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 100 -M /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}_markdups_metric_file.txt -I /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}_cleaned_sorted_rg.bam -O /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}_final.bam

rm /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}_cleaned_sorted_rg.bam

samtools index /lustre/scratch/jmanthey/01_certhia_genomics/01_bam_files/${SGE_TASK_ID}_final.bam
