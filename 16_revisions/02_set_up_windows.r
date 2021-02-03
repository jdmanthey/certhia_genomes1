options(scipen=999)
	project_directory <- "/lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/01_basic_filter/"
	directory_name <- "certhia_windows1"
	cluster <- "quanah"
	
	# read in reference index
	# filtered to only include genotyped chromosomes
	ref_index <- read.table("certhia_genome_filtered.fai", stringsAsFactors=F)
	
	# define window size
	window_size <- 100000
	
	# make directories
	dir.create(directory_name)
	
	# define intervals and write to helper files
	helpers <- c()
	for(a in 1:nrow(ref_index)) {
		
		a_start <- 1
		a_end <- a_start + window_size - 1
		a_max <- ref_index[a,2]
		a_windows <- ceiling((a_max - a_start) / window_size)
		a_chromosome <- ref_index[a,1]
		
		# loop for defining helper info for each window
		for(b in 1:a_windows) {
			if(b == a_windows) {
				a_end <- a_max
			}
			helpers <- rbind(helpers, c(paste(a_chromosome, ".recode.vcf.gz", sep=""), 
				paste(a_chromosome, ":", a_start, "-", a_end, sep="")))
			a_start <- a_start + window_size
			a_end <- a_end + window_size
		}
	}
	write(helpers[,1], file=paste(directory_name, "/helper9.txt", sep=""), ncolumns=1)
	write(helpers[,2], file=paste(directory_name, "/helper10.txt", sep=""), ncolumns=1)
	
	
	# write the array script
	a.script <- paste(directory_name, "/window_split_array.sh", sep="")
	write("#!/bin/sh", file=a.script)
	write("#SBATCH --chdir=./", file=a.script, append=T)
	write(paste("#SBATCH --job-name=", "window_split", sep=""), file=a.script, append=T)
	write("#SBATCH --nodes=1 --ntasks=1", file=a.script, append=T)
	write(paste("#SBATCH --partition ", cluster, sep=""), file=a.script, append=T)
	write("#SBATCH --time=02:00:00", file=a.script, append=T)
	write("#SBATCH --mem-per-cpu=8G", file=a.script, append=T)
	write(paste("#SBATCH --array=1-", nrow(helpers), sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("input_array=$( head -n${SLURM_ARRAY_TASK_ID} helper9.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("index_array=$( head -n${SLURM_ARRAY_TASK_ID} helper10.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	tabix_command <- paste("tabix ", project_directory, "${input_array} ${index_array} > ", 
		"temp${SLURM_ARRAY_TASK_ID}.vcf", sep="")
	write(tabix_command, file=a.script, append=T)
	write("", file=a.script, append=T)
	cat_command <- paste("cat ", project_directory, "header.vcf ", 
		"temp${SLURM_ARRAY_TASK_ID}.vcf > ", project_directory, "windows/${index_array}.vcf", sep="")
	write(cat_command, file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("rm temp${SLURM_ARRAY_TASK_ID}.vcf"), file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("bcftools query -f \'%POS\\t%REF\\t%ALT[\\t%GT]\\n\' ", 
		project_directory, "windows/${index_array}.vcf > ", 
		project_directory, "windows/${index_array}.simple.vcf", sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("rm ", project_directory, "windows/${index_array}.vcf", sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
