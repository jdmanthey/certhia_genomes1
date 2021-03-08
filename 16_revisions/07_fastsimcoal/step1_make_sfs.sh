
# repeat for file with no inversions
cd /lustre/scratch/jmanthey/05_certhia_genomics/04_fastsimcoal

# find total number of sites genotyped
grep -v '^#' fsc_no_inv.vcf | wc -l
# = 36265183

# remove invariant sites (file too big)
vcftools --vcf fsc_no_inv.vcf --min-alleles 2 --mac 1 --recode --recode-INFO-all --out fsc_no_inv

# find number of polymorphic sites genotyped
grep -v '^#' fsc_no_inv.recode.vcf | wc -l
# = 213997
# number of invariant sites = 36265183 - 213997 = 36051186

# convert vcf to sfs file
~/easySFS/easySFS.py -i fsc_no_inv.recode.vcf -p /lustre/scratch/jmanthey/05_certhia_genomics/04_fastsimcoal/sfs_popmap.txt -a -f --preview
~/easySFS/easySFS.py -i fsc_no_inv.recode.vcf -p /lustre/scratch/jmanthey/05_certhia_genomics/04_fastsimcoal/sfs_popmap.txt -a -f --proj 24,24

