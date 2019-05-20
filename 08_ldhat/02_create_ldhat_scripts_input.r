
# turn off scientific notation for output with large numbers in sliding windows
options(scipen=999)

snp_block <- 10000
chromosomes <- c("Ca_0002_Tg_1", "Ca_0005_Tg_1A", "Ca_0001_Tg_2", "Ca_0003_Tg_3", "Ca_0006_Tg_4", "Ca_0014_Tg_4A", "Ca_0007_Tg_5", "Ca_0009_Tg_6", "Ca_0008_Tg_7", "Ca_0010_Tg_8", "Ca_0011_Tg_9", "Ca_0015_Tg_10", "Ca_0012_Tg_11", "Ca_0013_Tg_12", "Ca_0016_Tg_13", "Ca_0019_Tg_14", "Ca_0020_Tg_15", "Ca_0022_Tg_17", "Ca_0021_Tg_18", "Ca_0025_Tg_19", "Ca_0029_Tg_19", "Ca_0017_Tg_20", "Ca_0026_Tg_21", "Ca_0028_Tg_23", "Ca_0023_Tg_24", "Ca_0004_Tg_Z")

dir.create("submit_scripts")
dir.create("input_files")


# loop for each chromosome scripts and input files
for(a in 1:length(chromosomes)) {
	
	# dump trash
	gc()
	
	# read files for snps on each chromosome
	a_temp2 <- read.table(paste(chromosomes[a], "_albescens.recode.vcf", sep=""))
	a_temp1 <- read.table(paste(chromosomes[a], "_americana.recode.vcf", sep=""))
		
	# how many jobs per this chromosome?
	n_jobs1 <- floor(nrow(a_temp1) / (snp_block - 200))
	n_jobs2 <- floor(nrow(a_temp2) / (snp_block - 200))
	if(n_jobs1 == 0) {
		n_jobs1 <- 1
	}
	if(n_jobs2 == 0) {
		n_jobs2 <- 1
	}
	
	# write job submission scripts
	base_name_1 <- paste("submit_scripts/americana_", chromosomes[a], ".sh", sep="")
	base_name_2 <- paste("submit_scripts/albescens_", chromosomes[a], ".sh", sep="")
	# write americana scripts
	write("#!/bin/sh", file=base_name_1, ncolumns=1)
	write("#$ -V", file=base_name_1, ncolumns=1, append=T)
	write("#$ -cwd", file=base_name_1, ncolumns=1, append=T)
	write("#$ -S /bin/bash", file=base_name_1, ncolumns=1, append=T)
	write(paste("#$ -N americana_", chromosomes[a], sep=""), file=base_name_1, ncolumns=1, append=T)
	write("#$ -q omni", file=base_name_1, ncolumns=1, append=T)
	write("#$ -pe sm 1", file=base_name_1, ncolumns=1, append=T)
	write("#$ -P quanah", file=base_name_1, ncolumns=1, append=T)
	write("#$ -l h_rt=48:00:00", file=base_name_1, ncolumns=1, append=T)
	write("#$ -l h_vmem=16G", file=base_name_1, ncolumns=1, append=T)
	write(paste("#$ -t 1-", n_jobs1, sep=""), file=base_name_1, ncolumns=1, append=T)
	write("", file=base_name_1, ncolumns=1, append=T)
	write(paste("/lustre/work/jmanthey/LDhat-master/interval -seq /lustre/scratch/jmanthey/ldhat/input_files/americana_", chromosomes[a], "_${SGE_TASK_ID}.sites -loc /lustre/scratch/jmanthey/ldhat/input_files/americana_", chromosomes[a], "_${SGE_TASK_ID}.locs  -lk /lustre/work/jmanthey/LDhat-master/input/certhia_lk.txt -bpen 10 -its 5000000 -samp 50000 -prefix am_", chromosomes[a], ".${SGE_TASK_ID}", sep=""), file=base_name_1, ncolumns=1, append=T)
	# write albescens scripts
	write("#!/bin/sh", file=base_name_2, ncolumns=1)
	write("#$ -V", file=base_name_2, ncolumns=1, append=T)
	write("#$ -cwd", file=base_name_2, ncolumns=1, append=T)
	write("#$ -S /bin/bash", file=base_name_2, ncolumns=1, append=T)
	write(paste("#$ -N albescens_", chromosomes[a], sep=""), file=base_name_2, ncolumns=1, append=T)
	write("#$ -q omni", file=base_name_2, ncolumns=1, append=T)
	write("#$ -pe sm 1", file=base_name_2, ncolumns=1, append=T)
	write("#$ -P quanah", file=base_name_2, ncolumns=1, append=T)
	write("#$ -l h_rt=48:00:00", file=base_name_2, ncolumns=1, append=T)
	write("#$ -l h_vmem=16G", file=base_name_2, ncolumns=1, append=T)
	write(paste("#$ -t 1-", n_jobs2, sep=""), file=base_name_2, ncolumns=1, append=T)
	write("", file=base_name_2, ncolumns=1, append=T)
	write(paste("/lustre/work/jmanthey/LDhat-master/interval -seq /lustre/scratch/jmanthey/ldhat/input_files/albescens_", chromosomes[a], "_${SGE_TASK_ID}.sites -loc /lustre/scratch/jmanthey/ldhat/input_files/albescens_", chromosomes[a], "_${SGE_TASK_ID}.locs  -lk /lustre/work/jmanthey/LDhat-master/input/certhia_lk.txt -bpen 10 -its 5000000 -samp 50000 -prefix al_", chromosomes[a], ".${SGE_TASK_ID}", sep=""), file=base_name_2, ncolumns=1, append=T)
		
	# genotypes lineage 1
	genotypes1 <- c()
	names1 <- c("Certhia_americana_KU128765", "Certhia_americana_KU128766", "Certhia_americana_KU128768",
	"Certhia_americana_KU128773", "Certhia_americana_KU128774", "Certhia_americana_KU128775", "Certhia_americana_KU29806",
	"Certhia_americana_KU29909", "Certhia_americana_KU31303", "Certhia_americana_UWBM113162", "Certhia_americana_UWBM113167", 
	"Certhia_americana_UWBM113168")
	for(b in 10:ncol(a_temp1)) {
		b_temp <- as.character(a_temp1[,b])
		b_temp <- sapply(strsplit(b_temp, ":"), "[[", 1)
		b_temp <- gsub("\\|", "/", b_temp)
		b_rep1 <- sapply(strsplit(b_temp, "/"), "[[", 1)
		b_rep2 <- sapply(strsplit(b_temp, "/"), "[[", 2)
		b_rep <- rep("N", length(b_rep1))
		b_rep[b_rep1 == "."] <- "?"
		b_rep[b_rep1 == "0" & b_rep1 == b_rep2] <- 0
		b_rep[b_rep1 == "1" & b_rep1 == b_rep2] <- 1
		b_rep[b_rep1 != b_rep2] <- 2
		genotypes1 <- c(genotypes1, paste(b_rep, collapse=""))
	}

	# genotypes lineage 2
	genotypes2 <- c()
	names2 <- c("Certhia_albescens_KU128978", "Certhia_albescens_UWBM106936", "Certhia_albescens_UWBM107890",
	 "Certhia_albescens_UWBM107897", "Certhia_albescens_KU128979", "Certhia_albescens_KU128986", 
	 "Certhia_albescens_KU131233", "Certhia_albescens_KU131234", "Certhia_albescens_KU131235", 
	 "Certhia_albescens_KU31300", "Certhia_albescens_KU31301", "Certhia_albescens_KU31309")
	for(b in 10:ncol(a_temp2)) {
		b_temp <- as.character(a_temp2[,b])
		b_temp <- sapply(strsplit(b_temp, ":"), "[[", 1)
		b_temp <- gsub("\\|", "/", b_temp)
		b_rep1 <- sapply(strsplit(b_temp, "/"), "[[", 1)
		b_rep2 <- sapply(strsplit(b_temp, "/"), "[[", 2)
		b_rep <- rep("N", length(b_rep1))
		b_rep[b_rep1 == "."] <- "?"
		b_rep[b_rep1 == "0" & b_rep1 == b_rep2] <- 0
		b_rep[b_rep1 == "1" & b_rep1 == b_rep2] <- 1
		b_rep[b_rep1 != b_rep2] <- 2
		genotypes2 <- c(genotypes2, paste(b_rep, collapse=""))
	}


	# loop for number of jobs in lineage 1
	start <- 1
	end <- snp_block
	for(b in 1:n_jobs1) {
		
		# base name of output files
		b_base_name <- paste("americana_", chromosomes[a], "_", b, sep="")
		
		genotypes <- c() # for .sites file
		distances <- c() # for .locs file
		locations <- c() # for mapping back to locations later
		
		# set up length of snps for this file
		if(b != n_jobs1) {
			n_snps <- snp_block
		} else { # length may be extra for last job for each chromosome
			n_snps <- nrow(a_temp1[start:nrow(a_temp1), ])
			end <- nrow(a_temp1)
		}
		
		# write sites file in lineage 1
		b_sites_header <- paste(ncol(a_temp1) - 9, n_snps, 2, sep="\t")
		write(b_sites_header, file=paste("input_files/", b_base_name, ".sites", sep=""), ncolumns=1)
		for(c in 1:length(names1)) {
			write(paste(">", names1[c], sep=""), file=paste("input_files/", b_base_name, ".sites", sep=""), ncolumns=1, append=T)
			write(substr(genotypes1[c], start, end), file=paste("input_files/", b_base_name, ".sites", sep=""), ncolumns=1, append=T)
		}
		
		# write locs file in lineage 1
		b_locs_header <- paste(n_snps, 1001, "L", sep="\t")
		b_locations1 <- a_temp1[start:end,2]
		b_total_length <- (b_locations1[length(b_locations1)] - b_locations1[1]) / 1000
		b_lengths <- (b_locations1 - b_locations1[1] + 1) / b_total_length
		write(b_locs_header, file=paste("input_files/", b_base_name, ".locs", sep=""), ncolumns=1)
		write(format(b_lengths, digits=10, scientific=F), file=paste("input_files/", b_base_name, ".locs", sep=""), ncolumns=1, append=T)
		
		# write mapping file (for remapping output of ldhat to specific locations)
		write(a_temp1[start:end, 2], file=paste("input_files/", b_base_name, ".mapping", sep=""), ncolumns=1)
		
		start <- start + snp_block - 200
		end <- end + snp_block - 200
	}
	
###################################	
###################################	
###################################	
###################################	
	# loop for number of jobs in lineage 2
	start <- 1
	end <- snp_block
	for(b in 1:n_jobs2) {
		
		# base name of output files
		b_base_name <- paste("albescens_", chromosomes[a], "_", b, sep="")
		
		genotypes <- c() # for .sites file
		distances <- c() # for .locs file
		locations <- c() # for mapping back to locations later
		
		# set up length of snps for this file
		if(b != n_jobs2) {
			n_snps <- snp_block
		} else { # length may be extra for last job for each chromosome
			n_snps <- nrow(a_temp2[start:nrow(a_temp2), ])
			end <- nrow(a_temp2)
		}
		
		# write sites file in lineage 1
		b_sites_header <- paste(ncol(a_temp2) - 9, n_snps, 2, sep="\t")
		write(b_sites_header, file=paste("input_files/", b_base_name, ".sites", sep=""), ncolumns=1)
		for(c in 1:length(names2)) {
			write(paste(">", names2[c], sep=""), file=paste("input_files/", b_base_name, ".sites", sep=""), ncolumns=1, append=T)
			write(substr(genotypes2[c], start, end), file=paste("input_files/", b_base_name, ".sites", sep=""), ncolumns=1, append=T)
		}
		
		# write locs file in lineage 1
		b_locs_header <- paste(n_snps, 1001, "L", sep="\t")
		b_locations2 <- a_temp2[start:end,2]
		b_total_length <- (b_locations2[length(b_locations2)] - b_locations2[1]) / 1000
		b_lengths <- (b_locations2 - b_locations2[1] + 1) / b_total_length
		write(b_locs_header, file=paste("input_files/", b_base_name, ".locs", sep=""), ncolumns=1)
		write(format(b_lengths, digits=10, scientific=F), file=paste("input_files/", b_base_name, ".locs", sep=""), ncolumns=1, append=T)
		
		# write mapping file (for remapping output of ldhat to specific locations)
		write(a_temp2[start:end, 2], file=paste("input_files/", b_base_name, ".mapping", sep=""), ncolumns=1)
		
		start <- start + snp_block - 200
		end <- end + snp_block - 200
	}

}


