# use a simplified vcf as input for creating a treemix file with snps at least minimum_dist apart

# column headers
# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	1	10	11	12	13	14	15	16	17	18	19	2	20	21	22	23	24	25	3	4	5	6	7	8	9
vcf_colnames <- c("POS", "REF", "ALT", "1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "2", "20", "21", "22", "23", "24", "25", "3", "4", "5", "6", "7", "8", "9")
options(scipen=999)

	# set up names of individuals
	ind_names <- c("C_alb_Chiricahua_KU128978", "C_alb_Mexico_UWBM106936", "C_alb_Mexico_UWBM107890", "C_alb_Mexico_UWBM107897", "C_ame_Pinaleno_KU128765", "C_ame_Pinaleno_KU128766", "C_ame_Pinaleno_KU128768", "C_ame_SCatalina_KU128773", "C_ame_SCatalina_KU128774", "C_ame_SCatalina_KU128775", "C_ame_Pinal_KU29806", "C_alb_Chiricahua_KU128979", "C_ame_Pinal_KU29909", "C_ame_Pinal_KU31303", "C_ame_Utah_UWBM113162", "C_ame_Utah_UWBM113167", "C_ame_Utah_UWBM113168", "C_fam_Outgroup_KU92846", "C_alb_Chiricahua_KU128986", "C_alb_Huachuca_KU131233", "C_alb_Huachuca_KU131234", "C_alb_Huachuca_KU131235", "C_alb_SRita_KU31300", "C_alb_SRita_KU31301", "C_alb_SRita_KU31309")

# minimum distance between snps
minimum_dist <- 10000

# set up population names
populations_labels <- sapply(strsplit(ind_names, "_"), "[[", 3)
populations_unique <- unique(populations_labels)
populations_unique <- populations_unique[c(6,5,3,4,1,9,8,2,7)] # order north to south

# write output 
write(populations_unique, file="certhia_contact_10000.treemix", sep=" ", ncolumns=9)

# all the files to read
x_files <- list.files(pattern="*simple.vcf")

# loop for each chromosome
for(a in 1:length(x_files)) {
	print(a)
	print(paste("Reading file"))
	a_rep <- read.table(x_files[a], stringsAsFactors=F, sep="\t")
	print(paste("Removing indels"))
	# remove indels
	a_rep <- a_rep[nchar(a_rep[,2]) == 1 & nchar(a_rep[,3]) == 1, ]
	print(paste("Removing variants with missing data"))
	# remove sites with missing data
	for(b in 4:ncol(a_rep)) {
		a_rep <- a_rep[grepl(pattern="\\./\\.", a_rep[,b]) == F,]
	}
	print(paste("Subsetting SNPs by distance"))
	# go through and decide which snps to keep
	new_a_rep <- a_rep[1,]
	a_rep <- a_rep[2:nrow(a_rep),]
	while(nrow(a_rep) > 1) {
		a_rep <- a_rep[abs(new_a_rep[nrow(new_a_rep), 1] - a_rep[,1]) >= minimum_dist, ]
		new_a_rep <- rbind(new_a_rep, a_rep[1,])
		a_rep <- a_rep[2:nrow(a_rep),]
	}
	if(nrow(a_rep) == 1) {
		new_a_rep <- rbind(new_a_rep, a_rep[1,])
	}
	new_a_rep <- na.omit(new_a_rep)

	print(paste("Writing output"))	
	# write output for each snp
	for(b in 1:nrow(new_a_rep)) {
		allele_freq_list <- c()
		b_rep <- as.character(new_a_rep[b,4:ncol(new_a_rep)])
		b_rep <- gsub("\\|", "/", b_rep)
		# loop for each population
		for(d in 1:length(populations_unique)) {
			# get this population's genotypes
			d_rep <- b_rep[populations_labels == populations_unique[d]]
			alt_allele <- length(d_rep[d_rep == "1/1"]) * 2 + length(d_rep[d_rep == "0/1"]) * 1
			ref_allele <- length(d_rep) * 2 - alt_allele
			allele_freq_list <- c(allele_freq_list, paste(ref_allele, alt_allele, sep=","))
		}
		write(allele_freq_list, file="certhia_contact_10000.treemix", sep=" ", ncolumns=9, append=T)
	}
}
