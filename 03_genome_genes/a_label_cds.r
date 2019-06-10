
fasta <- scan("certhia_cds.fasta", what="character")
gff <- read.table("certhia_round2_cds.gff", stringsAsFactors=F)

# extract the fasta file orientation for each sequence
fasta_orientation <- fasta[1:(length(fasta) / 2) *2 - 1]
fasta_orientation <- substr(fasta_orientation, nchar(fasta_orientation) - 2, nchar(fasta_orientation))

# extract the sequences
fasta_sequences <- fasta[1:(length(fasta) / 2) *2]

# gene names
gene_names <- sapply(strsplit(sapply(strsplit(gff[,9], "-RA:cds"), "[[", 1), "ID="), "[[", 2)
unique_gene_names <- unique(gene_names)

# loop for each gene
gene_sequences <- list()
orientations <- c()
for(a in 1:length(unique_gene_names)) {
	orientations <- c(orientations, fasta_orientation[gene_names == unique_gene_names[a]][1])
	gene_sequences[[a]] <- paste(fasta_sequences[gene_names == unique_gene_names[a]], collapse="")
}

# write the output
for(a in 1:length(unique_gene_names)) {
	a_name <- paste(">", unique_gene_names[a], sep="")
	if(a == 1) {
		write(a_name, file="certhia_cds_renamed.fasta", ncolumns=1)
	} else {
		write(a_name, file="certhia_cds_renamed.fasta", ncolumns=1, append=T)		
	}
	write(gene_sequences[[a]], file="certhia_cds_renamed.fasta", ncolumns=1, append=T)
}








# ficedula

fasta <- scan("ficedula_cds.fasta", what="character")
gff <- read.table("ficedula_cds.gff", stringsAsFactors=F, sep="\t", quote="\"")
fasta <- toupper(fasta)

# extract the sequences
fasta_sequences <- fasta[1:(length(fasta) / 2) *2]

# gene names
gene_names <- sapply(strsplit(sapply(strsplit(gff[,9], "gene="), "[[", 2), ";product"), "[[", 1)
unique_gene_names <- unique(gene_names)

# loop for each gene
gene_sequences <- list()
for(a in 1:length(unique_gene_names)) {
	gene_sequences[[a]] <- paste(fasta_sequences[gene_names == unique_gene_names[a]], collapse="")
}

# make the gene names the locus plus product
gene_names2 <- sapply(strsplit(sapply(strsplit(gff[match(unique_gene_names, gene_names),9], "product="), 
		"[[", 2), ";protein"), "[[", 1)
unique_gene_names <- paste(unique_gene_names, gene_names2, sep="_")
unique_gene_names <- gsub(';partial=true', '', unique_gene_names)

# write the output
for(a in 1:length(unique_gene_names)) {
	a_name <- paste(">", unique_gene_names[a], sep="")
	if(a == 1) {
		write(a_name, file="ficedula_cds_renamed.fasta", ncolumns=1)
	} else {
		write(a_name, file="ficedula_cds_renamed.fasta", ncolumns=1, append=T)		
	}
	write(gene_sequences[[a]], file="ficedula_cds_renamed.fasta", ncolumns=1, append=T)
}










# parus

fasta <- scan("parus_cds.fasta", what="character")
gff <- read.table("parus_cds.gff", stringsAsFactors=F, sep="\t", quote="\"")
fasta <- toupper(fasta)

# extract the sequences
fasta_sequences <- fasta[1:(length(fasta) / 2) *2]

# gene names
gene_names <- sapply(strsplit(sapply(strsplit(gff[,9], "gene="), "[[", 2), ";product"), "[[", 1)
unique_gene_names <- unique(gene_names)

# loop for each gene
gene_sequences <- list()
for(a in 1:length(unique_gene_names)) {
	gene_sequences[[a]] <- paste(fasta_sequences[gene_names == unique_gene_names[a]], collapse="")
}

# make the gene names the locus plus product
gene_names2 <- sapply(strsplit(sapply(strsplit(gff[match(unique_gene_names, gene_names),9], "product="), 
		"[[", 2), ";protein"), "[[", 1)
unique_gene_names <- paste(unique_gene_names, gene_names2, sep="_")
unique_gene_names <- gsub(';partial=true', '', unique_gene_names)

# write the output
for(a in 1:length(unique_gene_names)) {
	a_name <- paste(">", unique_gene_names[a], sep="")
	if(a == 1) {
		write(a_name, file="parus_cds_renamed.fasta", ncolumns=1)
	} else {
		write(a_name, file="parus_cds_renamed.fasta", ncolumns=1, append=T)		
	}
	write(gene_sequences[[a]], file="parus_cds_renamed.fasta", ncolumns=1, append=T)
}









# taeniopygia

fasta <- scan("taeniopygia_cds.fasta", what="character")
gff <- read.table("taeniopygia_cds.gff", stringsAsFactors=F, sep="\t", quote="\"")
fasta <- toupper(fasta)

# extract the sequences
fasta_sequences <- fasta[1:(length(fasta) / 2) *2]

# gene names
gene_names <- sapply(strsplit(sapply(strsplit(gff[,9], "gene="), "[[", 2), ";product"), "[[", 1)
unique_gene_names <- unique(gene_names)

# loop for each gene
gene_sequences <- list()
for(a in 1:length(unique_gene_names)) {
	gene_sequences[[a]] <- paste(fasta_sequences[gene_names == unique_gene_names[a]], collapse="")
}

# make the gene names the locus plus product
gene_names2 <- sapply(strsplit(sapply(strsplit(gff[match(unique_gene_names, gene_names),9], "product="), 
		"[[", 2), ";protein"), "[[", 1)
unique_gene_names <- paste(unique_gene_names, gene_names2, sep="_")
unique_gene_names <- gsub(';partial=true', '', unique_gene_names)

# write the output
for(a in 1:length(unique_gene_names)) {
	a_name <- paste(">", unique_gene_names[a], sep="")
	if(a == 1) {
		write(a_name, file="taeniopygia_cds_renamed.fasta", ncolumns=1)
	} else {
		write(a_name, file="taeniopygia_cds_renamed.fasta", ncolumns=1, append=T)		
	}
	write(gene_sequences[[a]], file="taeniopygia_cds_renamed.fasta", ncolumns=1, append=T)
}








