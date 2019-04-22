require(Biostrings)
genome <- readDNAStringSet("05_certhia_hic.fasta")
genome2 <- genome
genome_names <- genome@ranges@NAMES

# synteny matches 
x <- read.table("certhia_filtered.txt", sep="\t", stringsAsFactors=F, header=T)
scaffolds_matched <- unique(x[,11])
finch_chromosomes <- c("Tg_1", "Tg_1A", "Tg_2", "Tg_3", "Tg_4", "Tg_4A", "Tg_5", "Tg_6", "Tg_7", "Tg_8", "Tg_9", "Tg_10", "Tg_11", "Tg_12", "Tg_13", "Tg_14", "Tg_15", "Tg_17", "Tg_18", "Tg_19", "Tg_20", "Tg_21", "Tg_22", "Tg_23", "Tg_24", "Tg_26", "Tg_27", "Tg_28", "Tg_Z")         

# setup opposite of in function
"%ni%" <- Negate("%in%")

# find matches vs zebra finch genome for each scaffold
sorting_names <- c()
for(a in 1:length(scaffolds_matched)) {
	a_rep <- x[x[,11] == scaffolds_matched[a],]
	a_matches <- sort(table(a_rep[,10]), decreasing=T)
	# if greater than 70% matches to a zebra finch chromosome, match it with that, otherwise = "none"
	if(as.numeric(a_matches[1]) / nrow(a_rep) > 0.7) {
		a_match <- names(a_matches)[1]
	} else {
		a_match <- "none"
	}
	sorting_names <- c(sorting_names, a_match)
}
finch_chromosomes <- finch_chromosomes[finch_chromosomes %in% sorting_names] 
no_matches <- scaffolds_matched[sorting_names == "none"]

# make a new genome file that is renamed and reordered

# first add the reordered scaffolds for big zebra finch chromosomes
for(a in 1:length(finch_chromosomes)) {
	# matching scaffold name
	a_original <- scaffolds_matched[sorting_names == finch_chromosomes[a]]
	# select that scaffold
	a_genome <- genome2[genome2@ranges@NAMES %in% a_original]
	# remove the one scaffold from the overall object
	genome2 <- genome2[genome2@ranges@NAMES %ni% a_original]
	# rename the scaffold
	a_original <- substr(a_original, 14, nchar(a_original))
	for(b in 1:length(a_original)) {
		if(nchar(a_original[b]) == 1) {
			a_original[b] <- paste("000", a_original[b], sep="")
		} else if(nchar(a_original[b]) == 2) {
			a_original[b] <- paste("00", a_original[b], sep="")
		} else if(nchar(a_original[b]) == 3) {
			a_original[b] <- paste("0", a_original[b], sep="")
		}
	}
	a_name <- paste("Ca_", a_original, "_", finch_chromosomes[a], sep="")
	a_genome@ranges@NAMES <- a_name
	if(a == 1) {
		new_genome <- a_genome
	} else {
		new_genome <- c(new_genome, a_genome)
	}
}

# second add the big scaffolds that didn't clearly match a big zebra finch scaffold
for(a in 1:length(no_matches)) {
	a_original <- no_matches[a]
	# select that scaffold
	a_genome <- genome2[genome2@ranges@NAMES %in% a_original]
	# remove the one scaffold from the overall object
	genome2 <- genome2[genome2@ranges@NAMES %ni% a_original]
	# rename the scaffold
	a_original <- substr(a_original, 14, nchar(a_original))
	if(nchar(a_original) == 1) {
		a_original <- paste("000", a_original, sep="")
	} else if(nchar(a_original) == 2) {
		a_original[b] <- paste("00", a_original, sep="")
	} else if(nchar(a_original) == 3) {
		a_original <- paste("0", a_original, sep="")
	}
	a_name <- paste("Ca_", a_original, sep="")
	a_genome@ranges@NAMES <- a_name
	new_genome <- c(new_genome, a_genome)
}

# now add all of the smaller scaffolds that also need to be renamed
a_original <- genome2@ranges@NAMES
a_original <- substr(a_original, 14, nchar(a_original))
for(a in 1:length(a_original)) {
	if(nchar(a_original[a]) == 1) {
		a_original[a] <- paste("000", a_original[a], sep="")
	} else if(nchar(a_original[a]) == 2) {
		a_original[a] <- paste("00", a_original[a], sep="")
	} else if(nchar(a_original[a]) == 3) {
		a_original[a] <- paste("0", a_original[a], sep="")
	}
}
a_name <- paste("Ca_", a_original, sep="")
genome2@ranges@NAMES <- a_name
new_genome <- c(new_genome, genome2)

# write output
writeXStringSet(new_genome, file="06_certhia_reordered.fasta")



