#!/bin/sh
#SBATCH --chdir=./
#SBATCH --job-name=filter
#SBATCH --nodes=1 --ntasks=1
#SBATCH --partition quanah
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=1-26

# define input files from helper file during genotyping
input_array=$( head -n${SLURM_ARRAY_TASK_ID} vcf_list.txt | tail -n1 )
input_array=${input_array%.g.vcf}

# define main working directory
workdir=/lustre/scratch/jmanthey/05_certhia_genomics/03_vcf

# run vcftools with basic filtering
vcftools --vcf ${workdir}/${input_array}.g.vcf --max-missing 0.6 --minQ 20 --minGQ 20 --minDP 5 --max-meanDP 50 --remove-indels --recode --recode-INFO-all --out ${workdir}/01_basic_filter/${input_array}
# run bcftools to simplify the vcftools output
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/01_basic_filter/${input_array}.recode.vcf > ${workdir}/01_basic_filter/${input_array}.simple.vcf

# run vcftools with basic filtering and maf of 0.05
vcftools --vcf ${workdir}/${input_array}.g.vcf --max-missing 0.6 --minQ 20 --minGQ 20 --minDP 5 --max-meanDP 50 --maf 0.05 --remove-indels --recode --recode-INFO-all --out ${workdir}/02_basic_filter_maf/${input_array}
# run bcftools to simplify the vcftools output
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/02_basic_filter_maf/${input_array}.recode.vcf > ${workdir}/02_basic_filter_maf/${input_array}.simple.vcf

# run vcftools for biallelic SNPs separated by 10kbp with no missing data
vcftools --vcf ${workdir}/${input_array}.g.vcf --max-missing 1.0 --minQ 20 --minGQ 20 --minDP 5 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 1 --thin 10000 --remove-indels --recode --recode-INFO-all --out ${workdir}/03_biallelic_10kbp/${input_array}
# run bcftools to simplify the vcftools output
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/03_biallelic_10kbp/${input_array}.recode.vcf > ${workdir}/03_biallelic_10kbp/${input_array}.simple.vcf

# run vcftools for biallelic SNPs separated by 50kbp with no missing data
vcftools --vcf ${workdir}/${input_array}.g.vcf --max-missing 1.0 --minQ 20 --minGQ 20 --minDP 5 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 1 --thin 50000 --remove-indels --recode --recode-INFO-all --out ${workdir}/04_biallelic_50kbp/${input_array}
# run bcftools to simplify the vcftools output
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/04_biallelic_50kbp/${input_array}.recode.vcf > ${workdir}/04_biallelic_50kbp/${input_array}.simple.vcf

# run vcftools with basic filtering and only keeping biallelic SNPs
vcftools --vcf ${workdir}/${input_array}.g.vcf --max-missing 0.6 --minQ 20 --minGQ 20 --minDP 5 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 1 --remove-indels --recode --recode-INFO-all --out ${workdir}/05_biallelic_40p/${input_array}
# run bcftools to simplify the vcftools output
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/05_biallelic_40p/${input_array}.recode.vcf > ${workdir}/05_biallelic_40p/${input_array}.simple.vcf

# run vcftools for biallelic SNPs separated by 50kbp with no missing data and maf of 0.2
vcftools --vcf ${workdir}/${input_array}.g.vcf --max-missing 1.0 --minQ 20 --minGQ 20 --minDP 5 --max-meanDP 50 --min-alleles 2 --max-alleles 2 --mac 1 --maf 0.2 --thin 50000 --remove-indels --recode --recode-INFO-all --out ${workdir}/06_LD/${input_array}
# run bcftools to simplify the vcftools output
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n ' ${workdir}/06_LD/${input_array}.recode.vcf > ${workdir}/06_LD/${input_array}.simple.vcf
