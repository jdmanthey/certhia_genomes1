# make sure output directory written already 
# in hpcc script: mkdir output

abba_certhia <- function(xxx, output_number) {
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

	# set up names and numbers for the taxa in these tests	
	abba_tests <- list()
	abba_tests[[1]] <- c("SCatalina", "Utah", "Mexico", "Outgroup")
	abba_tests[[2]] <- c("Pinaleno", "Utah", "Mexico", "Outgroup")
	abba_tests[[3]] <- c("Pinal", "Utah", "Mexico", "Outgroup")
	abba_tests[[4]] <- c("Chiricahua", "Mexico", "Utah", "Outgroup")
	abba_tests[[5]] <- c("Huachuca", "Mexico", "Utah", "Outgroup")
	abba_tests[[6]] <- c("SRita", "Mexico", "Utah", "Outgroup")
			
	# write initial output file
	write(paste("pop1", "pop2", "pop3", "n_snps", "d", "z-score", sep="\t"), file=paste("output/abba_", output_number, ".txt", sep=""), ncolumns=6)
	
	# subsets for each set of populations and remove sites with missing data
	# and keep only sites biallelic in the quartet of pops	
	for(a in 1:length(abba_tests)) {
		# the four populations subsets of the data frame
		a_rep1 <- x[,group_numbers[group_names == abba_tests[[a]][1]]]
		a_rep2 <- x[,group_numbers[group_names == abba_tests[[a]][2]]]
		a_rep3 <- x[,group_numbers[group_names == abba_tests[[a]][3]]]
		a_rep4 <- x[,group_numbers[group_names == abba_tests[[a]][4]]]
		
		# group all together to check if these loci are biallelic in the designated populations
		a_rep_total <- cbind(a_rep1, a_rep2, a_rep3, a_rep4)
		biallelic <- apply(a_rep_total, 1, find_biallelic)
		
		# subset to biallelic
		a_rep1 <- a_rep1[biallelic, ]
		a_rep2 <- a_rep2[biallelic, ]
		a_rep3 <- a_rep3[biallelic, ]
		a_rep4 <- a_rep4[biallelic]
		a_rep_total <- a_rep_total[biallelic, ]
		
		# from the grouped sample, find loci that have data for each population
		non_missing <- apply(a_rep_total, 1, find_missing)
		
		# subset the snps without missing data
		a_rep1 <- a_rep1[non_missing, ]
		a_rep2 <- a_rep2[non_missing, ]
		a_rep3 <- a_rep3[non_missing, ]
		a_rep4 <- a_rep4[non_missing]
		a_rep_total <- a_rep_total[non_missing, ]
		
		# find proportions derived allele for each snp and population
		d1 <- c()
		d2 <- c()
		d3 <- c()
		# loop for each snp (coulnd't figure out how to use apply with a matrix and vector to match the outgroup)
		for(b in 1:nrow(a_rep1)) {
			b_rep1 <- as.character(a_rep1[b,])
			b_rep2 <- as.character(a_rep2[b,])
			b_rep3 <- as.character(a_rep3[b,])
			
			# remove missing data
			b_rep1 <- b_rep1[b_rep1 != "./."]
			b_rep2 <- b_rep2[b_rep2 != "./."]
			b_rep3 <- b_rep3[b_rep3 != "./."]
			
			# proportion each allele
			b_rep1 <- (length(b_rep1[b_rep1 == a_rep4[b]]) * 2 + length(b_rep1[b_rep1 == "0/1"])) / (length(b_rep1) * 2)
			b_rep2 <- (length(b_rep2[b_rep2 == a_rep4[b]]) * 2 + length(b_rep2[b_rep2 == "0/1"])) / (length(b_rep2) * 2)
			b_rep3 <- (length(b_rep3[b_rep3 == a_rep4[b]]) * 2 + length(b_rep3[b_rep3 == "0/1"])) / (length(b_rep3) * 2)
			
			# add to output
			d1 <- c(d1, b_rep1)
			d2 <- c(d2, b_rep2)
			d3 <- c(d3, b_rep3)
		}
		
		# d calculation
		a_rep_D <- sum( (1 - d1) * d2 * d3 * 1 - d1 * (1 - d2) * d3 * 1 ) / sum( (1 - d1) * d2 * d3 * 1 + d1 * (1 - d2) * d3 * 1 )
		
		# 100 randomizations
		a_rep_D_random <- c()
		for(b in 1:100) {
			b_rep_random <- t(apply(a_rep_total[,1:9], 1, function(xxxx) {return(sample(xxxx, 9))}))
			# loop for each snp
			d1 <- c()
			d2 <- c()
			d3 <- c()
			for(c in 1:nrow(b_rep_random)) {
				b_rep1 <- as.character(b_rep_random[c,1:3])
				b_rep2 <- as.character(b_rep_random[c,4:6])
				b_rep3 <- as.character(b_rep_random[c,7:9])
			
				# remove missing data
				b_rep1 <- b_rep1[b_rep1 != "./."]
				b_rep2 <- b_rep2[b_rep2 != "./."]
				b_rep3 <- b_rep3[b_rep3 != "./."]
			
				# proportion each allele
				b_rep1 <- (length(b_rep1[b_rep1 == a_rep4[b]]) * 2 + length(b_rep1[b_rep1 == "0/1"])) / (length(b_rep1) * 2)
				b_rep2 <- (length(b_rep2[b_rep2 == a_rep4[b]]) * 2 + length(b_rep2[b_rep2 == "0/1"])) / (length(b_rep2) * 2)
				b_rep3 <- (length(b_rep3[b_rep3 == a_rep4[b]]) * 2 + length(b_rep3[b_rep3 == "0/1"])) / (length(b_rep3) * 2)
			
				# add to output
				d1 <- c(d1, b_rep1)
				d2 <- c(d2, b_rep2)
				d3 <- c(d3, b_rep3)
			}
			
			# d calculation
			a_rep_D_random <- c(a_rep_D_random, 
				sum( (1 - d1) * d2 * d3 * 1 - d1 * (1 - d2) * d3 * 1 ) / sum( (1 - d1) * d2 * d3 * 1 + d1 * (1 - d2) * d3 * 1 ))
		}
		
		# remove randomizations that didn't work because of only sampling some pops. with missing data
		#a_rep_D_random <- as.numeric(na.omit(a_rep_D_random))[1:100]
		
		# calculate z score
		z_rep <- (a_rep_D - 0) / sd(a_rep_D_random)
		
		# combine and then write the output
		output_rep <- c(abba_tests[[a]][1], abba_tests[[a]][2], abba_tests[[a]][3], nrow(a_rep_total), a_rep_D, z_rep)
		write(output_rep, file=paste("output/abba_", output_number, ".txt", sep=""), sep="\t", ncolumns=6, append=T)
	}
}




# find bi-allelic function
# use with apply over rows 
find_biallelic <- function(xxxx) {
	x_rep <- xxxx[xxxx != "./."]
	return(length(unique(x_rep)) > 1)
}

# find missingfunction
# use with apply over rows 
# pop 1, 2, 3 each have 3 individuals and pop4 has one individual
# for removal of any snp with no individuals in a population
find_missing <- function(xxxx) {
	x_rep1 <- xxxx[1:3]
	x_rep2 <- xxxx[4:6]
	x_rep3 <- xxxx[7:9]
	x_rep4 <- xxxx[10]
	x_rep1 <- x_rep1[x_rep1 != "./."]
	x_rep2 <- x_rep2[x_rep2 != "./."]
	x_rep3 <- x_rep3[x_rep3 != "./."]
	x_rep4 <- x_rep4[x_rep4 != "./."]
	x_missing_total <- length(x_rep1) + length(x_rep2) + length(x_rep3)
	if(length(x_rep1) == 0 | length(x_rep2) == 0 | length(x_rep3) == 0 | length(x_rep4) == 0 | x_missing_total < 7) {
		return(FALSE)
	} else {
		return(TRUE)
	}
}









