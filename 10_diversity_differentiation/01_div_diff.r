# make sure output directory written already 
# in hpcc script: mkdir output

div_certhia <- function(xxx, output_number) {
	# read in vcf window
	x <- read.table(xxx, sep="\t", stringsAsFactors=F)
	
	vcf_colnames <- c("POS", "REF", "ALT", "1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "2", "20", "21", "22", "23", "24", "25", "3", "4", "5", "6", "7", "8", "9")

	options(scipen=999)
	
	# set up names of individuals
	ind_names <- c("C_alb_Chiricahua_KU128978", "C_alb_Mexico_UWBM106936", "C_alb_Mexico_UWBM107890", "C_alb_Mexico_UWBM107897", "C_ame_Pinaleno_KU128765", "C_ame_Pinaleno_KU128766", "C_ame_Pinaleno_KU128768", "C_ame_SCatalina_KU128773", "C_ame_SCatalina_KU128774", "C_ame_SCatalina_KU128775", "C_ame_Pinal_KU29806", "C_alb_Chiricahua_KU128979", "C_ame_Pinal_KU29909", "C_ame_Pinal_KU31303", "C_ame_Utah_UWBM113162", "C_ame_Utah_UWBM113167", "C_ame_Utah_UWBM113168", "C_fam_Outgroup_KU92846", "C_alb_Chiricahua_KU128986", "C_alb_Huachuca_KU131233", "C_alb_Huachuca_KU131234", "C_alb_Huachuca_KU131235", "C_alb_SRita_KU31300", "C_alb_SRita_KU31301", "C_alb_SRita_KU31309")
	
	# set up population numbers
	group_names <- sapply(strsplit(ind_names, "_"), "[[", 3)
	# columns in vcf of these names
	group_numbers <- 1:length(group_names) + 3
	
	# set up different groups for differentiation tests
	parental_north <- group_numbers[group_names == "Utah"]
	parental_south <- group_numbers[group_names == "Mexico"]
	contact_north <- group_numbers[group_names == "SCatalina" | group_names == "Pinaleno" | group_names == "Pinal"]
	contact_south <- group_numbers[group_names == "Chiricahua" | group_names == "Huachuca" | group_names == "SRita"]
	outgroup <- group_numbers[group_names == "Outgroup"]
	total_north <- group_numbers[group_names == "SCatalina" | group_names == "Pinaleno" | group_names == "Pinal" | group_names == "Utah"]
	total_south <- group_numbers[group_names == "Chiricahua" | group_names == "Huachuca" | group_names == "SRita" | group_names == "Mexico"]
	
	# write initial output file
	write(c("test", "group", "n_sites", "n_variants", "stat"), file=paste("output/div_certhia_", output_number, ".txt", sep=""), ncolumns=5, sep="\t")
	
	# diversity for each individual
	for(a in 1:length(ind_names)) {
		# subset for this individual (skip first three lines)
		a_rep <- x[,a + 3]
		# remove missing genotypes
		a_rep <- a_rep[a_rep != "./."]
		# number of genotyped sites and number heterozygous
		a_sites <- length(a_rep)
		a_het_sites <- length(a_rep[a_rep == "0/1"])
		a_het <- a_het_sites / a_sites
		
		# combine and write diversity stats to output
		a_output <- c("heterozygosity", ind_names[a], a_sites, a_het_sites, a_het)
		write(a_output, file=paste("output/div_certhia_", output_number, ".txt", sep=""), ncolumns=5, sep="\t", append=T)
	}
	
	# set up dxy tests
	dxy_tests1 <- list()
	dxy_tests2 <- list()
	dxy_test_names <- list()
	dxy_tests1[[1]] <- parental_north
	dxy_tests2[[1]] <- parental_south
	dxy_test_names[[1]] <- "parentals_north_south" 
	dxy_tests1[[2]] <- contact_north
	dxy_tests2[[2]] <- contact_south
	dxy_test_names[[2]] <- "contact_north_south" 
	dxy_tests1[[3]] <- total_north
	dxy_tests2[[3]] <- total_south
	dxy_test_names[[3]] <- "total_north_south" 
	dxy_tests1[[4]] <- parental_north
	dxy_tests2[[4]] <- outgroup
	dxy_test_names[[4]] <- "parental_north_outgroup" 
	dxy_tests1[[5]] <- parental_south
	dxy_tests2[[5]] <- outgroup
	dxy_test_names[[5]] <- "parental_south_outgroup" 
	# loop for each dxy test
	for(a in 1:length(dxy_tests1)) {
		# subset data
		a_rep1 <- x[,dxy_tests1[[a]]]
		a_rep2 <- x[,dxy_tests2[[a]]]
		if(is.null(nrow(a_rep2))) {
			a_rep2 <- as.matrix(a_rep2)
		}
		
		# remove all sites with less than 3 individuals genotyped (or 1 in the case of the outgroup)
		a_keep1 <- apply(a_rep1, 1, filter_alleles)
		a_keep2 <- apply(a_rep2, 1, filter_alleles)
			#filtering
		a_rep1 <- a_rep1[a_keep1 & a_keep2, ]
		a_rep2 <- a_rep2[a_keep1 & a_keep2, ]
		if(is.null(nrow(a_rep2))) {
			a_rep2 <- as.matrix(a_rep2)
		}
		
		# count the total number of included genotyped sites at this point
		a_sites <- nrow(a_rep1)
		
		# combine the two matrices to find sites that are variant between the two pops
		a_total <- cbind(a_rep1, a_rep2)
		# find and keep variant sites
		a_variants <- apply(a_total, 1, find_variants)
		a_rep1 <- a_rep1[a_variants, ]
		a_rep2 <- a_rep2[a_variants, ]
		if(is.null(nrow(a_rep2))) {
			a_rep2 <- as.matrix(a_rep2)
		}
		
		# count the number of sites
		a_variant_sites <- nrow(a_rep1)
		
		# loop for each sites to calculate dxy
		dxy_all <- c()
		for(b in 1:nrow(a_rep1)) {
			# subset to this snp
			b_rep1 <- as.character(a_rep1[b,])
			b_rep2 <- as.character(a_rep2[b,])
			
			# measure proportion of reference allele 
			b_ref1 <- (length(b_rep1[b_rep1 == "0/0"]) * 2 + length(b_rep1[b_rep1 == "0/1"]) * 1) / (length(b_rep1) * 2)
			b_ref2 <- (length(b_rep2[b_rep2 == "0/0"]) * 2 + length(b_rep2[b_rep2 == "0/1"]) * 1) / (length(b_rep2) * 2)
			
			# calc dxy
			dxy_rep <- b_ref1 * (1 - b_ref2) + b_ref2 * (1 - b_ref1)
			dxy_all <- c(dxy_all, dxy_rep)

		}
		
		# sum of all dxy divided by the total number of variant sites
		a_dxy <- sum(dxy_all) / a_variant_sites
		
		# combine and write dxy stats to output
		a_output <- c("dxy", dxy_test_names[[a]], a_sites, a_variant_sites, a_dxy)
		write(a_output, file=paste("output/div_certhia_", output_number, ".txt", sep=""), ncolumns=5, sep="\t", append=T)
			
	}
	
	
	
	
	# set up fst tests
	fst_tests1 <- list()
	fst_tests2 <- list()
	fst_test_names <- list()
	fst_tests1[[1]] <- parental_north
	fst_tests2[[1]] <- parental_south
	fst_test_names[[1]] <- "parentals_north_south" 
	fst_tests1[[2]] <- contact_north
	fst_tests2[[2]] <- contact_south
	fst_test_names[[2]] <- "contact_north_south" 
	fst_tests1[[3]] <- total_north
	fst_tests2[[3]] <- total_south
	fst_test_names[[3]] <- "total_north_south" 

	# loop for each fst test
	for(a in 1:length(fst_tests1)) {
		# subset data
		a_rep1 <- x[,fst_tests1[[a]]]
		a_rep2 <- x[,fst_tests2[[a]]]
		
		# remove all sites with less than 3 individuals genotyped (or 1 in the case of the outgroup)
		a_keep1 <- apply(a_rep1, 1, filter_alleles)
		a_keep2 <- apply(a_rep2, 1, filter_alleles)
			#filtering
		a_rep1 <- a_rep1[a_keep1 & a_keep2, ]
		a_rep2 <- a_rep2[a_keep1 & a_keep2, ]
				
		# count the total number of included genotyped sites at this point
		a_sites <- nrow(a_rep1)
		
		# combine the two matrices to find sites that are variant between the two pops
		a_total <- cbind(a_rep1, a_rep2)
		# find and keep variant sites
		a_variants <- apply(a_total, 1, find_variants)
		a_rep1 <- a_rep1[a_variants, ]
		a_rep2 <- a_rep2[a_variants, ]
		if(is.null(nrow(a_rep2))) {
			a_rep2 <- as.matrix(a_rep2)
		}
		
		# count the number of sites
		a_variant_sites <- nrow(a_rep1)
		
		# loop for each site to calculate fst
		numerator_fst_all <- c()
		denominator_fst_all <- c()
		for(b in 1:nrow(a_rep1)) {
			# subset to this snp
			b_rep1 <- as.character(a_rep1[b,])
			b_rep2 <- as.character(a_rep2[b,])
			
			# fst is the reich et al. 2009 estimator for small sample sizes
			# equation presented nicer in Willing et al. 2012 page 9
			pop1_ind_count <- length(b_rep1) 
			pop2_ind_count <- length(b_rep2)
			alt_allele_count1 <- (2 * length(b_rep1[b_rep1 == "1/1"]) + 1 * length(b_rep1[b_rep1 == "0/1"]))
			alt_allele_count2 <- (2 * length(b_rep2[b_rep2 == "1/1"]) + 1 * length(b_rep2[b_rep2 == "0/1"]))
			all_allele_count1 <- 2 * length(b_rep1)
			all_allele_count2 <- 2 * length(b_rep2)
			expected_het1 <- (alt_allele_count1 * (all_allele_count1 - alt_allele_count1)) / 
			(all_allele_count1 * (all_allele_count1 - 1))
			expected_het2 <- (alt_allele_count2 * (all_allele_count2 - alt_allele_count2)) / 
			(all_allele_count2 * (all_allele_count2 - 1))
			
			# find the fst numerator and denominator values for this snp (they all get summed and divided for 
			# the final estimate)
			numerator_rep <- (alt_allele_count1 / (2 * pop1_ind_count) - 
			alt_allele_count2 / (2 * pop2_ind_count))^2 - (expected_het1 / (2 * pop1_ind_count)) - 
			(expected_het2 / (2 * pop2_ind_count))
			
			denominator_rep <- numerator_rep + expected_het1 + expected_het2
			
			# add to total outputs
			numerator_fst_all <- c(numerator_fst_all, numerator_rep)
			denominator_fst_all <- c(denominator_fst_all, denominator_rep)
		}
		
		# mean of all fst 
		a_fst <- sum(numerator_fst_all) / sum(denominator_fst_all) 
		
		# combine and write fst stats to output
		a_output <- c("fst", fst_test_names[[a]], a_sites, a_variant_sites, a_fst)
		write(a_output, file=paste("output/div_certhia_", output_number, ".txt", sep=""), ncolumns=5, sep="\t", append=T)
			
	}	
	
	
	
	
}


# find sites with variants
# function to use with apply across rows
find_variants <- function(xxxx) {
	xxxx <- xxxx[xxxx != "./."]
	xxxx <- unique(xxxx)
	if(length(xxxx) == 1) {
		if(xxxx != "0/1") {
			return(FALSE)
		} else {
			return(TRUE)
		}
	} else {
		return(TRUE)
	}
}


#filter alleles function
# across all rows using apply
# filtered for at least 3 individuals with the exception of the outgroup
# which is allowed no missing data since it is one individual
filter_alleles <- function(xxxx) {
	if(length(xxxx) == 1) {
		xxxx <- xxxx[xxxx != "./."]
		if(length(xxxx) == 1) {
			return(TRUE)
		} else {
			return(FALSE)
		}
	} else {
		xxxx <- xxxx[xxxx != "./."]
		if(length(xxxx) >= 3) {
			return(TRUE)
		} else {
			return(FALSE)
		}
	}
}



















