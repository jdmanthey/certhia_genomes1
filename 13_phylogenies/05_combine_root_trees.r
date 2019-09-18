require(ape)

x_files <- list.files(pattern="*tre$")

x_files2 <- sapply(strsplit(x_files, ".tre"), "[[", 1)
position <- list()
chromosome <- list()
# get the chromosome and position for each file
for(a in 1:length(x_files2)) {
	if(length(strsplit(x_files2[a], "_")[[1]]) == 5) {
		position[[a]] <- strsplit(x_files2[a], "_")[[1]][5]
		chromosome[[a]] <- paste(strsplit(x_files2[a], "_")[[1]][1], "_", strsplit(x_files2[a], "_")[[1]][2], "_", strsplit(x_files2[a], "_")[[1]][3], "_",
			strsplit(x_files2[a], "_")[[1]][4], sep="")
	} else if(length(strsplit(x_files2[a], "_")[[1]]) == 3) {
		position[[a]] <- strsplit(x_files2[a], "_")[[1]][3]
		chromosome[[a]] <- paste(strsplit(x_files2[a], "_")[[1]][1], "_", strsplit(x_files2[a], "_")[[1]][2], sep="")
	} 
}
x_locations <- unlist(position)
x_chrom <- unlist(chromosome)


# read in and root the trees
for(a in 1:length(x_files)) {
	a_rep <- read.tree(file=x_files[a])
	a_rep <- root(a_rep, outgroup="C_familiaris_KU92846", resolve.root=T)
	if(a == 1) {
		output_trees <- a_rep
	} else {
		output_trees <- c(output_trees, a_rep)
	}
}

# write the trees to a new file
write.tree(output_trees, file="certhia_total.trees")
write.table(cbind(seq(from=1, to=length(x_files), by=1), x_chrom, x_locations), file="certhia_total_notes.txt", quote=F, row.names=F, col.names=F)

