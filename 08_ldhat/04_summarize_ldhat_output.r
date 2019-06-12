# each snp block is 2000, with 200 overlapping the previous file
# need to remove the first 100 bp of most files and last 100 bp of most files

# this summarizes the output as estimates of 4Nr per kbp

# read in original mapping file and results file
maps_files <- list.files("input_files", pattern="*mapping", full.names=T)
results_files <- list.files("raw_results", pattern="*res.txt", full.names=T)

# info about the mapping and results files
mapping_chromosome <- paste(sapply(strsplit(maps_files, "_"), "[[", 3), sapply(strsplit(maps_files, "_"), "[[", 4), sapply(strsplit(maps_files, "_"), "[[", 5), sapply(strsplit(maps_files, "_"), "[[", 6), sep="_")
mapping_species <- substr(sapply(strsplit(sapply(strsplit(maps_files, "_"), "[[", 2), "/"), "[[", 2), 1, 2)
mapping_replicate <- sapply(strsplit(sapply(strsplit(maps_files, "\\.mapping"), "[[", 1), "_"), "[[", 7)

results_chromosome <- sapply(strsplit(paste(sapply(strsplit(results_files, "_"), "[[", 3), sapply(strsplit(results_files, "_"), "[[", 4), sapply(strsplit(results_files, "_"), "[[", 5), sapply(strsplit(results_files, "_"), "[[", 6), sep="_"), "\\."), "[[", 1)
results_species <- sapply(strsplit(sapply(strsplit(results_files, "_"), "[[", 2), "/"), "[[", 2)
results_replicate <- sapply(strsplit(sapply(strsplit(results_files, "\\."), "[[", 2), "res"), "[[", 1)



# americana first
# loop for each chromosome
write(paste("#ldhat recombination estimates"), file="ldhat_output.txt", sep="\t")
write(paste("species\tchromosome\tlocation\tmean_rho"), file="ldhat_output.txt", sep="\t", append=T)
for(a in 1:length(unique(results_chromosome))) {
	a_reps <- results_replicate[results_species == "am" & results_chromosome == unique(results_chromosome)[a]]
	output <- c()
	for(b in 1:length(a_reps)) {
		b_map <- read.table(paste("input_files/americana_", unique(results_chromosome)[a], "_", b, ".mapping", sep=""), sep="\t")
		b_results <- read.table(paste("raw_results/am_", unique(results_chromosome)[a], ".", b, "res.txt", sep=""), sep="\t", header=T)
		# divide by the size of the whole window (because ldhat assumes the units are kbp, but they are not)
		b_results[,2] <- b_results[,2] / ((b_map[nrow(b_map),1] - b_map[1,1] + 1) / 1000)
		b_output <- cbind(rep("am", nrow(b_map)), rep(unique(results_chromosome)[a], nrow(b_map)), b_map[,1], b_results[,2])
		if(b == 1 & length(a_reps) == 1) {
			# remove first base estimate of recombination (highly inflated) and keep rest of file
			b_output <- b_output[2:nrow(b_output),]
		} else if (b == 1) {
			# remove first base estimate of recombination (highly inflated) and remove final 100 bp of overlap
			b_output <- b_output[2:(nrow(b_output) - 100),]
		} else if (b == length(a_reps)) {
			# remove first 100 bp overlap and keep rest of file
			b_output <- b_output[101:nrow(b_output),]
		} else {
			# remove first and last 100 bp of overlap
			b_output <- b_output[101:(nrow(b_output) - 100),]
		}
		write.table(b_output, file="ldhat_output.txt", sep="\t", row.names=F, col.names=F, quote=F, append=T)
	}
}


# albescens second
for(a in 1:length(unique(results_chromosome))) {
	a_reps <- results_replicate[results_species == "al" & results_chromosome == unique(results_chromosome)[a]]
	output <- c()
	for(b in 1:length(a_reps)) {
		b_map <- read.table(paste("input_files/albescens_", unique(results_chromosome)[a], "_", b, ".mapping", sep=""), sep="\t")
		b_results <- read.table(paste("raw_results/al_", unique(results_chromosome)[a], ".", b, "res.txt", sep=""), sep="\t", header=T)
		# divide by the size of the whole window (because ldhat assumes the units are kbp, but they are not)
		b_results[,2] <- b_results[,2] / ((b_map[nrow(b_map),1] - b_map[1,1] + 1) / 1000)
		b_output <- cbind(rep("al", nrow(b_map)), rep(unique(results_chromosome)[a], nrow(b_map)), b_map[,1], b_results[,2])
		if(b == 1 & length(a_reps) == 1) {
			# remove first base estimate of recombination (highly inflated) and keep rest of file
			b_output <- b_output[2:nrow(b_output),]
		} else if (b == 1) {
			# remove first base estimate of recombination (highly inflated) and remove final 100 bp of overlap
			b_output <- b_output[2:(nrow(b_output) - 100),]
		} else if (b == length(a_reps)) {
			# remove first 100 bp overlap and keep rest of file
			b_output <- b_output[101:nrow(b_output),]
		} else {
			# remove first and last 100 bp of overlap
			b_output <- b_output[101:(nrow(b_output) - 100),]
		}
		write.table(b_output, file="ldhat_output.txt", sep="\t", row.names=F, col.names=F, quote=F, append=T)
	}
}


# summarize into windows
window_size <- 100000

options(scipen=999)

x <- read.table("ldhat_output.txt", sep="\t", header=T, stringsAsFactors=F)
# chromosome order
chromosomes <- c("Ca_0002_Tg_1", "Ca_0005_Tg_1A", "Ca_0001_Tg_2", "Ca_0003_Tg_3", "Ca_0006_Tg_4", "Ca_0014_Tg_4A", "Ca_0007_Tg_5", "Ca_0009_Tg_6", "Ca_0008_Tg_7", "Ca_0010_Tg_8", "Ca_0011_Tg_9", "Ca_0015_Tg_10", "Ca_0012_Tg_11", "Ca_0013_Tg_12", "Ca_0016_Tg_13", "Ca_0019_Tg_14", "Ca_0020_Tg_15", "Ca_0022_Tg_17", "Ca_0021_Tg_18", "Ca_0025_Tg_19", "Ca_0029_Tg_19", "Ca_0017_Tg_20", "Ca_0026_Tg_21", "Ca_0028_Tg_23", "Ca_0023_Tg_24", "Ca_0004_Tg_Z")

# write output
write(c("chr", "start", "end", "window", "americana_n_snps", "americana_rho", "albescens_n_snps", "albescens_rho"), file="certhia_rho_windows.txt", sep="\t", ncolumns=8)

# window counter
window_counter <- 1
# loop for each chromosome
for(a in 1:length(chromosomes)) {
	# subset data
	a_rep <- x[x$chromosome == chromosomes[a],]
	# number of windows
	a_num_windows <- floor(max(a_rep$location) / window_size)
	
	a_start <- 1
	a_end <- window_size
	# loop for each window
	for(b in 1:a_num_windows) {
		b_albescens <- a_rep[a_rep$species == "al" & a_rep$location >= a_start & a_rep$location <= a_end, ]
		b_americana <- a_rep[a_rep$species == "am" & a_rep$location >= a_start & a_rep$location <= a_end, ]
		b_output <- c(chromosomes[a], a_start, a_end, window_counter, nrow(b_americana), 
			mean(b_americana$mean_rho), nrow(b_albescens), mean(b_albescens$mean_rho))
		write(b_output, file="certhia_rho_windows.txt", sep="\t", ncolumns=8, append=T)
		window_counter <- window_counter + 1
		a_start <- a_start + window_size
		a_end <- a_end + window_size
	}
	
	# reduce size of entire matrix
	x <- x[x$chromosome != chromosomes[a],]
}
































