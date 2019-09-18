

options(scipen=999)

# define the number of chromosomes (i.e., the number of sets of smc files)
chr_number <- 25

# number of bootstraps 
bootstraps <- 10

# length of chunks to sample
chunks <- 100000

# create output directories
for(a in 1:bootstraps) {
	dir.create(paste("bootstrap_", a, sep=""))
}

# loop for each chromosome
for(a in 1:chr_number) {
	# read in files
	a_rep1 <- read.table(paste(a, "_a_total.smc", sep=""), sep=" ", stringsAsFactors=F)
	a_rep2 <- read.table(paste(a, "_b_total.smc", sep=""), sep=" ", stringsAsFactors=F)
	a_rep3 <- read.table(paste(a, "_a_americana.smc", sep=""), sep=" ", stringsAsFactors=F)
	a_rep4 <- read.table(paste(a, "_b_americana.smc", sep=""), sep=" ", stringsAsFactors=F)
	a_rep5 <- read.table(paste(a, "_a_albescens.smc", sep=""), sep=" ", stringsAsFactors=F)
	a_rep6 <- read.table(paste(a, "_b_albescens.smc", sep=""), sep=" ", stringsAsFactors=F)
	
	# get cumulative sums of the positions for each of the files
	cs1 <- cumsum(a_rep1[,1])
	cs2 <- cumsum(a_rep2[,1])
	cs3 <- cumsum(a_rep3[,1])
	cs4 <- cumsum(a_rep4[,1])
	cs5 <- cumsum(a_rep5[,1])
	cs6 <- cumsum(a_rep6[,1])
	
	# total length of chromosome
	max_length <- max(cumsum(a_rep2[,1]))
	
	# number of chunks to sample 
	number_chunks <- floor(max_length / chunks)
	
	# loop for each bootstrap
	for(b in 1:bootstraps) {
		# sample random chunks
		start_of_sampling_chunks <- sample(seq(from=1, to=(max_length - chunks), by=1), size=number_chunks, replace=T)
		end_of_sampling_chunks <- start_of_sampling_chunks + chunks
		
		# loop for each set of sampling chunks
		for(d in 1:length(start_of_sampling_chunks)) {
			# write the output
			if(d == 1) {
				# file 1
				write.table(a_rep1[cs1 > start_of_sampling_chunks[d] & cs1 < end_of_sampling_chunks[d],], 
					file=paste("bootstrap_", b, "/", a, "_a_total.smc", sep=""), sep=" ", quote=F, row.names=F, col.names=F)
				write.table(a_rep2[cs2 > start_of_sampling_chunks[d] & cs2 < end_of_sampling_chunks[d],], 
					file=paste("bootstrap_", b, "/", a, "_b_total.smc", sep=""), sep=" ", quote=F, row.names=F, col.names=F)
				write.table(a_rep3[cs3 > start_of_sampling_chunks[d] & cs3 < end_of_sampling_chunks[d],], 
					file=paste("bootstrap_", b, "/", a, "_a_americana.smc", sep=""), sep=" ", quote=F, row.names=F, col.names=F)			
				write.table(a_rep4[cs4 > start_of_sampling_chunks[d] & cs4 < end_of_sampling_chunks[d],], 
					file=paste("bootstrap_", b, "/", a, "_b_americana.smc", sep=""), sep=" ", quote=F, row.names=F, col.names=F)	
				write.table(a_rep5[cs5 > start_of_sampling_chunks[d] & cs5 < end_of_sampling_chunks[d],], 
					file=paste("bootstrap_", b, "/", a, "_a_albescens.smc", sep=""), sep=" ", quote=F, row.names=F, col.names=F)
				write.table(a_rep6[cs6 > start_of_sampling_chunks[d] & cs6 < end_of_sampling_chunks[d],], 
					file=paste("bootstrap_", b, "/", a, "_b_albescens.smc", sep=""), sep=" ", quote=F, row.names=F, col.names=F)
			} else {
				write.table(a_rep1[cs1 > start_of_sampling_chunks[d] & cs1 < end_of_sampling_chunks[d],], 
					file=paste("bootstrap_", b, "/", a, "_a_total.smc", sep=""), sep=" ", quote=F, row.names=F, col.names=F, append=T)
				write.table(a_rep2[cs2 > start_of_sampling_chunks[d] & cs2 < end_of_sampling_chunks[d],], 
					file=paste("bootstrap_", b, "/", a, "_b_total.smc", sep=""), sep=" ", quote=F, row.names=F, col.names=F, append=T)
				write.table(a_rep3[cs3 > start_of_sampling_chunks[d] & cs3 < end_of_sampling_chunks[d],], 
					file=paste("bootstrap_", b, "/", a, "_a_americana.smc", sep=""), sep=" ", quote=F, row.names=F, col.names=F, append=T)			
				write.table(a_rep4[cs4 > start_of_sampling_chunks[d] & cs4 < end_of_sampling_chunks[d],], 
					file=paste("bootstrap_", b, "/", a, "_b_americana.smc", sep=""), sep=" ", quote=F, row.names=F, col.names=F, append=T)	
				write.table(a_rep5[cs5 > start_of_sampling_chunks[d] & cs5 < end_of_sampling_chunks[d],], 
					file=paste("bootstrap_", b, "/", a, "_a_albescens.smc", sep=""), sep=" ", quote=F, row.names=F, col.names=F, append=T)
				write.table(a_rep6[cs6 > start_of_sampling_chunks[d] & cs6 < end_of_sampling_chunks[d],], 
					file=paste("bootstrap_", b, "/", a, "_b_albescens.smc", sep=""), sep=" ", quote=F, row.names=F, col.names=F, append=T)				
			}
		}
	}
}
