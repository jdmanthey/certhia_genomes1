
options(scipen=999)

x_files <- list.files(pattern="*sites")
x_names <- paste(substr(x_files, 1, nchar(x_files) - 6), "_mask.bed", sep="")

# loop for each mask
for(a in 1:length(x_files)) {
	# read in sites
	a_rep <- read.table(x_files[a], sep="\t")
	# get chromosome
	a_chrom <- as.character(a_rep[1,1])
	# get sites and append zero
	a_sites <- c(0,as.numeric(a_rep[,2]))
	# get difference between sites
	a_diff <- diff(a_sites)
	# create two indexes for the a_diff object 
	a_index1 <- seq(from=1, to=length(a_diff), by=1) - 1
	a_index2 <- seq(from=1, to=length(a_diff), by=1)
	# remove the zero from the a_sites_object
	a_sites <- a_sites[2:length(a_sites)]
	
	# subset the indexes by the a_diff object 
	a_index1 <- a_index1[a_diff != 1]
	if(a_index1[1] == 0) {a_index1[1] <- 1}
	a_index2 <- a_index2[a_diff != 1]
	
	
	# keep info on sites that are more than 1bp apart (two objects = beginning of gap and end of gap)
	a_sites1 <- a_sites[a_index1]
	a_sites2 <- a_sites[a_index2]
	if(a_sites1[1] == a_sites2[1]) { a_sites1[1] <- 0 }
	# the a_sites1 object is already in the right format, because it includes that last sequenced base (0-based)
	# need to adjust the second base to the position that is before the next sequenced base in the gap (1-based)
	a_sites2 <- a_sites2 - 1
	
	#make the output bed
	output <- cbind(rep(a_chrom, length(a_sites2)), a_sites1, a_sites2)
	write.table(output, file=x_names[a], quote=F, sep="\t", row.names=F, col.names=F)
}
