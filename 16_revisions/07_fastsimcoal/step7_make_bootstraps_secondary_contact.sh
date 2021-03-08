# determine number of genotype lines in files
grep -v "^#" fsc_no_inv.recode.vcf | wc -l
# 213997 in file w/o inversions

# split the files into 500 SNP blocks (and also grep out headers)
grep "^#" fsc_no_inv.recode.vcf > header.sites
grep -v "^#" fsc_no_inv.recode.vcf > fsc_no_inv.sites
mkdir split
split -l 500 fsc_no_inv.sites split/x

# generate 100 new files with 500 SNP bootstrap blocks
cd 02_no_inv/
for i in {1..100}; do
	# define output name
	out_name=fsc_no_inv__${i}.vcf
	# sfs name
	sfs_name=fsc_no_inv__${i}_jointMAFpop1_0.obs
	# new sfs name
	new_sfs_name=bs${i}/secondary_contact_jointMAFpop1_0.obs
	# make directory
	mkdir bs${i}
	
	# get all blocks of SNPs
	# add 428 blocks of SNPs to get close to total number of SNPs (214000) 
	for r in {1..428}; do cat `find ./split -type f -name "*" | shuf -n1` >> temp.txt; done
	# change location for each so that none are duplicates (easySFS removes those)
	awk '$2=NR' OFS="\t" temp.txt > temp2.txt
	
	# add the header
	cat header.sites > $out_name
	# add the data
	cat temp2.txt >> $out_name
	
	# run easysfs
	~/easySFS/easySFS.py -i $out_name -p /lustre/scratch/jmanthey/05_certhia_genomics/04_fastsimcoal/sfs_popmap.txt -a -f --proj 24,24

	# move the sfs file to the working directory
	mv output/fastsimcoal2/$sfs_name ./$new_sfs_name

	# remove the output sfs directory and the intermediate vcf files
	rm -r output
	rm $out_name
	rm temp.txt
	rm temp2.txt

	# Progress
 	echo bs$i" done"
done



# add number of invariant sites genotyped to each SFS
for i in {1..100}; do 
	sed -i '/^d1_0/s/\t0/\t36051186/1' bs${i}/secondary_contact_jointMAFpop1_0.obs
done


