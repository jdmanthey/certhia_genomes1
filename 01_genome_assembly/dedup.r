x_files <- list.files(pattern="split_merge*")

for(a in 1:length(x_files)) {
	a_rep <- read.table(x_files[a], sep=" ", stringsAsFactors=F)
	# remove identical alignments
	a_rep1 <- paste(a_rep[,2], a_rep[,3], a_rep[,6], a_rep[,7])
	a_rep <- a_rep[match(unique(a_rep1), a_rep1),]
	
	# map quality >= 10
	a_rep <- a_rep[a_rep[,9] >= 10 & a_rep[,12] >= 10, ]
	# remove lines that have less than a 500 bp alignment distance
	a_rep1 <- a_rep[(a_rep[,2] == a_rep[,6] & abs(a_rep[,3] - a_rep[,7]) < 500) == FALSE, ]
	
	# create output file
	write.table(a_rep1, file=paste("no_dupes_", a, ".txt", sep=""), sep=" ", quote=F, row.names=F, col.names=F)


}



x_files <- list.files(pattern="no_dupes_*")

# read in files
input <- list()
for(a in 1:length(x_files)) {
	print(a)
	input[[a]] <- read.table(paste("no_dupes_", a, ".txt", sep=""), sep=" ", stringsAsFactors=F)
}
# total number of alignments
sum(unlist(lapply(input, nrow)))

# combine all the matrices
input2 <- do.call(rbind, input)

# remove identical alignments
a_rep1 <- paste(input2[,2], input2[,3], input2[,6], input2[,7])
input2 <- input2[match(unique(a_rep1), a_rep1),]

# check mapping quality for both reads (columns 9 and 12)
input2 <- input2[input2[,9] >= 10 & input2[,12] >= 10, ]

# write output of deduped total file
write.table(input2, file="merged_nodups.txt", sep=" ", quote=F, row.names=F, col.names=F)

