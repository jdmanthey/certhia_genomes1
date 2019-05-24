# subset the vcf files into sliding windows and keep a record of which file is which

# column headers
# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	1	10	11	12	13	14	15	16	17	18	19	2	20	21	22	23	24	25	3	4	5	6	7	8	9
vcf_colnames <- c("POS", "REF", "ALT", "1", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "2", "20", "21", "22", "23", "24", "25", "3", "4", "5", "6", "7", "8", "9")
options(scipen=999)

	# set up names of individuals
	ind_names <- c("C_alb_Chiricahua_KU128978", "C_alb_Mexico_UWBM106936", "C_alb_Mexico_UWBM107890", "C_alb_Mexico_UWBM107897", "C_ame_Pinaleno_KU128765", "C_ame_Pinaleno_KU128766", "C_ame_Pinaleno_KU128768", "C_ame_SCatalina_KU128773", "C_ame_SCatalina_KU128774", "C_ame_SCatalina_KU128775", "C_ame_Pinal_KU29806", "C_alb_Chiricahua_KU128979", "C_ame_Pinal_KU29909", "C_ame_Pinal_KU31303", "C_ame_Utah_UWBM113162", "C_ame_Utah_UWBM113167", "C_ame_Utah_UWBM113168", "C_fam_Outgroup_KU92846", "C_alb_Chiricahua_KU128986", "C_alb_Huachuca_KU131233", "C_alb_Huachuca_KU131234", "C_alb_Huachuca_KU131235", "C_alb_SRita_KU31300", "C_alb_SRita_KU31301", "C_alb_SRita_KU31309")
	# set up population numbers
	group_names <- sapply(strsplit(ind_names, "_"), "[[", 3)
	# columns in vcf of these names
	group_numbers <- 1:length(group_names) + 3
	
file_number <- 1

# list the files and get the chromosome names
x_files <- list.files(pattern="*.simple.vcf")
# order the files by chr number in the zebra finch
x_files <- x_files[c(2,5,1,3,6,14,7,9,8,10,11,15,12,13,16,19,20,22,21,25,29,17,26,28,23,4)]
x_names <- sapply(strsplit(x_files, ".simple"), "[[", 1)

# make window directory
dir.create("window_vcf")

# define window size
window_size <- 100000

# define minimum number of sites genotyped per window to do the tests
min_variants <- 30000

# output tracking file
output_tracker <- c()

for(a in 3:length(x_files)) {
	print(a)
	print(paste("Reading Table"))
	a_start <- 1
	a_end <- window_size
	a_rep <- read.table(x_files[a], stringsAsFactors=F, sep="\t")
	colnames(a_rep) <- vcf_colnames
	
	print(paste("Filtering"))
	# remove indels and sites with more than two alleles
	a_rep <- a_rep[nchar(a_rep$REF) == 1 & nchar(a_rep$ALT) == 1, ]

	
	print(paste("Modification"))	
	# replace all phased haplotypes with the regular /
	for(b in 4:ncol(a_rep)) {
		a_rep[,b] <- gsub("\\|", "/", a_rep[,b])
	}
	
	# define the number of windows
	a_windows <- floor(a_rep[nrow(a_rep), 1] / window_size)
	
	print(paste("Looping windows"))
	# loop for each window
	for(b in 1:a_windows) {
		b_rep <- a_rep[a_rep[,1] >= a_start & a_rep[,1] <= a_end, ]
		if(nrow(b_rep) >= min_variants) {
			write.table(b_rep, file=paste("window_vcf/", file_number, ".txt", sep=""), quote=F, sep="\t", row.names=F, col.names=F)
			# update the output tracking file
			b_output <- c(x_names[a], a_start, a_end, file_number)
			output_tracker <- rbind(output_tracker, b_output)
		}
		
		# move the sliding window
		a_start <- a_start + window_size
		a_end <- a_end + window_size
		# update the file number
		file_number <- file_number + 1
	}
	gc()
}

write.table(output_tracker, file="div_diff_windows_info.txt", sep="\t", quote=F, row.names=F, col.names=F)


# write r scripts for submitting jobs
for(a in 1:file_number) {
	a_name <- paste("window_vcf/div_diff_", a, ".r", sep="")
	a_rep <- 'source("01_div_diff.r")'
	write(a_rep, a_name, ncolumns=1)
	a_rep <- paste('div_certhia("', a, '.txt", ', a, ")", sep="")
	write(a_rep, a_name, ncolumns=1, append=T)
}





