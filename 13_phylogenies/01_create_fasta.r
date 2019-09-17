
# genotype all scaffolds > 2 Mbp
chr_list_to_genotype <- c("Ca_0002_Tg_1", "Ca_0005_Tg_1A", "Ca_0001_Tg_2", "Ca_0003_Tg_3", "Ca_0006_Tg_4", "Ca_0014_Tg_4A", "Ca_0007_Tg_5", "Ca_0009_Tg_6", "Ca_0008_Tg_7", "Ca_0010_Tg_8", "Ca_0011_Tg_9", "Ca_0015_Tg_10", "Ca_0012_Tg_11", "Ca_0013_Tg_12", "Ca_0016_Tg_13", "Ca_0019_Tg_14", "Ca_0020_Tg_15", "Ca_0022_Tg_17", "Ca_0021_Tg_18", "Ca_0025_Tg_19", "Ca_0029_Tg_19", "Ca_0017_Tg_20", "Ca_0026_Tg_21", "Ca_0028_Tg_23", "Ca_0023_Tg_24", "Ca_0004_Tg_Z", "Ca_0018", "Ca_0024", "Ca_0027", "Ca_0030")

x_files <- list.files(pattern="*phylo.simple.vcf")
x_names <- sapply(strsplit(x_files, "_phylo.simple.vcf"), "[[", 1)

# window size
window_size <- 50000
# minimum content to keep the window (e.g., 0.5 = half the window size)
keep_window <- 0.5
# minimum window content per individual to keep the window 
ind_window <- 0.1

output_directory <- "certhia_fasta"
dir.create(output_directory)

# order of individuals in vcf: 
# 1	10	11	12	13	14	15	16	17	18	19	2	20	21	22	23	24	25	3	4	5	6	7	8	9

# set up names of individuals
ind_names <- c("C_albescens_KU128978", "C_albescens_UWBM106936", "C_albescens_UWBM107890", "C_albescens_UWBM107897",
	"C_americana_KU128765", "C_americana_KU128766", "C_americana_KU128768", "C_americana_KU128773", "C_americana_KU128774",
	"C_americana_KU128775", "C_americana_KU29806", "C_albescens_KU128979", "C_americana_KU29909", "C_americana_KU31303",
	"C_americana_UWBM113162", "C_americana_UWBM113167", "C_americana_UWBM113168", "C_familiaris_KU92846", "C_albescens_KU128986", 
	"C_albescens_KU131233", "C_albescens_KU131234", "C_albescens_KU131235", "C_albescens_KU31300", "C_albescens_KU31301", 
	"C_albescens_KU31309")
ind_fasta_names <- paste(">", ind_names, sep="")


# loop for each chromosome
for(a in 1:length(x_files)) {
	a_start <- 1
	a_end <- window_size
	a_rep <- read.table(x_files[a], sep="\t", stringsAsFactors=F)
	# remove non bi-allelic SNPs (how to code so many possibilities and null alleles?)
	a_rep <- a_rep[nchar(a_rep[,3]) == 1, ]
	a_window_number <- ceiling(max(a_rep[,1]) / window_size)
	
	# loop for each window
	for(b in 1:a_window_number) {
		b_rep <- a_rep[a_rep[,1] >= a_start & a_rep[,1] <= a_end, ]
		
		# check missing sites per individual
		keep_going <- TRUE
		for(d in 4:ncol(b_rep)) {
			d_test <- b_rep[,d]
			d_test <- d_test[d_test != "./."]
			if(length(d_test) < (window_size * ind_window)) {
				print(length(d_test))
				keep_going <- FALSE
			}
		}
		
		if(nrow(b_rep) > (window_size * keep_window) & keep_going == TRUE) {
			allele1 <- b_rep[,2]
			allele2 <- b_rep[,3]
			# set up heterozygous codes
			heterozygous <- rep("?", length(allele1))
			heterozygous[allele1 == "A" & allele2 == "C"] <- "M"
			heterozygous[allele1 == "C" & allele2 == "A"] <- "M"
			heterozygous[allele1 == "A" & allele2 == "G"] <- "R"
			heterozygous[allele1 == "G" & allele2 == "A"] <- "R"
			heterozygous[allele1 == "A" & allele2 == "T"] <- "W"
			heterozygous[allele1 == "T" & allele2 == "A"] <- "W"
			heterozygous[allele1 == "C" & allele2 == "G"] <- "S"
			heterozygous[allele1 == "G" & allele2 == "C"] <- "S"
			heterozygous[allele1 == "C" & allele2 == "T"] <- "Y"
			heterozygous[allele1 == "T" & allele2 == "C"] <- "Y"
			heterozygous[allele1 == "G" & allele2 == "T"] <- "K"
			heterozygous[allele1 == "T" & allele2 == "G"] <- "K"
				
			# loop for each individual
			for(c in 1:length(ind_names)) {
				c_rep <- b_rep[,c+3]
				c_rep[c_rep == "./."] <- "?"
				c_rep[c_rep == "0/0"] <- allele1[c_rep == "0/0"]
				c_rep[c_rep == "0|0"] <- allele1[c_rep == "0|0"]
				c_rep[c_rep == "1/1"] <- allele2[c_rep == "1/1"]
				c_rep[c_rep == "1|1"] <- allele2[c_rep == "1|1"]
				c_rep[c_rep == "0/1"] <- heterozygous[c_rep == "0/1"]
				c_rep[c_rep == "0|1"] <- heterozygous[c_rep == "0|1"]
				c_rep[c_rep == "1|0"] <- heterozygous[c_rep == "1|0"]
				c_rep <- paste(c_rep, collapse="", sep="")
				if(c == 1) {
					write(ind_fasta_names[c], paste(output_directory, "/", x_names[a], "_", a_start, ".fasta", sep=""), ncolumns=1)
					write(c_rep, paste(output_directory, "/", x_names[a], "_", a_start, ".fasta", sep=""), ncolumns=1, append=T)
				} else {
					write(ind_fasta_names[c], paste(output_directory, "/", x_names[a], "_", a_start, ".fasta", sep=""), ncolumns=1, append=T)
					write(c_rep, paste(output_directory, "/", x_names[a], "_", a_start, ".fasta", sep=""), ncolumns=1, append=T)
				}
			}
		}
		a_start <- a_start + window_size
		a_end <- a_end + window_size
	}
	
	
	
}





