

tabix /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/Ca_0002_Tg_1.recode.vcf.gz \
Ca_0002_Tg_1:45600001-64600000 > temp1.vcf

cat /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/header.vcf \
temp1.vcf > /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0002_Tg_1:45600001-64600000.vcf

rm temp1.vcf

vcftools --vcf /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0002_Tg_1:45600001-64600000.vcf \
--max-missing 1.0 --maf 0.04 --recode --recode-INFO-all --remove-indv 25 --min-alleles 2 --max-alleles 2 --mac 1 \
--out /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0002_Tg_1:45600001-64600000

bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' \
/lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0002_Tg_1:45600001-64600000.recode.vcf > \
/lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0002_Tg_1:45600001-64600000.simple.vcf





tabix /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/Ca_0006_Tg_4.recode.vcf.gz \
Ca_0006_Tg_4:800001-5600000 > temp1.vcf

cat /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/header.vcf \
temp1.vcf > /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0006_Tg_4:800001-5600000.vcf

rm temp1.vcf

vcftools --vcf /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0006_Tg_4:800001-5600000.vcf \
--max-missing 1.0 --maf 0.04 --recode --recode-INFO-all --remove-indv 25 --min-alleles 2 --max-alleles 2 --mac 1 \
--out /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0006_Tg_4:800001-5600000

bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' \
/lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0006_Tg_4:800001-5600000.recode.vcf > \
/lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0006_Tg_4:800001-5600000.simple.vcf





tabix /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/Ca_0008_Tg_7.recode.vcf.gz \
Ca_0008_Tg_7:28100001-36000000 > temp1.vcf

cat /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/header.vcf \
temp1.vcf > /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0008_Tg_7:28100001-36000000.vcf

rm temp1.vcf

vcftools --vcf /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0008_Tg_7:28100001-36000000.vcf \
--max-missing 1.0 --maf 0.04 --recode --recode-INFO-all --remove-indv 25 --min-alleles 2 --max-alleles 2 --mac 1 \
--out /lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0008_Tg_7:28100001-36000000

bcftools query -f '%POS\t%REF\t%ALT[\t%GT]\n' \
/lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0008_Tg_7:28100001-36000000.recode.vcf > \
/lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/05_biallelic_40p/outliers/Ca_0008_Tg_7:28100001-36000000.simple.vcf


