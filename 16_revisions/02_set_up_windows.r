	options(scipen=999)
	project_directory <- "/lustre/scratch/jmanthey/05_certhia_genomics/03_vcf/02_basic_filter_maf/"
	directory_name <- "certhia_windows_maf"
	cluster <- "quanah"
	max_number_jobs <- 50
	
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
	
	# calculate number of array jobs
	if(nrow(helpers) > max_number_jobs) {
		n_array_jobs <- max_number_jobs
		n_jobs_per_array <- ceiling(nrow(helpers) / max_number_jobs)
	} else {
		n_array_jobs <- nrow(helpers)
		n_jobs_per_array <- 1
	}

	# write the array script
	a.script <- paste(directory_name, "/window_split_array.sh", sep="")
	write("#!/bin/sh", file=a.script)
	write("#SBATCH --chdir=./", file=a.script, append=T)
	write(paste("#SBATCH --job-name=", "window_split", sep=""), file=a.script, append=T)
	write("#SBATCH --nodes=1 --ntasks=1", file=a.script, append=T)
	write(paste("#SBATCH --partition ", cluster, sep=""), file=a.script, append=T)
	write("#SBATCH --time=48:00:00", file=a.script, append=T)
	write("#SBATCH --mem-per-cpu=8G", file=a.script, append=T)
	write(paste("#SBATCH --array=1-", n_array_jobs, sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	
	write("# Set the number of runs that each SLURM task should do", file=a.script, append=T)
	write(paste("PER_TASK=", n_jobs_per_array, sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	
	write("# Calculate the starting and ending values for this task based", file=a.script, append=T)
	write("# on the SLURM task and the number of runs per task.", file=a.script, append=T)
	write("START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))", file=a.script, append=T)
	write("END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))", file=a.script, append=T)
	write("", file=a.script, append=T)
	
	write("# Print the task and run range", file=a.script, append=T)
	write("echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM", file=a.script, append=T)
	write("", file=a.script, append=T)

	write("# Run the loop of runs for this task.", file=a.script, append=T)	
	write("for (( run=$START_NUM; run<=END_NUM; run++ )); do", file=a.script, append=T)
	write("\techo This is SLURM task $SLURM_ARRAY_TASK_ID, run number $run", file=a.script, append=T)
	write("", file=a.script, append=T)
	
	write("\tinput_array=$( head -n${run} helper9.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	write("\tindex_array=$( head -n${run} helper10.txt | tail -n1 )", file=a.script, append=T)
	write("", file=a.script, append=T)
	tabix_command <- paste("\ttabix ", project_directory, "${input_array} ${index_array} > ", 
		"temp${run}.vcf", sep="")
	write(tabix_command, file=a.script, append=T)
	write("", file=a.script, append=T)
	cat_command <- paste("\tcat ", project_directory, "header.vcf ", 
		"temp${run}.vcf > ", project_directory, "windows/${index_array}.vcf", sep="")
	write(cat_command, file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("\trm temp${run}.vcf"), file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("\tbcftools query -f \'%POS\\t%REF\\t%ALT[\\t%GT]\\n\' ", 
		project_directory, "windows/${index_array}.vcf > ", 
		project_directory, "windows/${index_array}.simple.vcf", sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write(paste("\trm ", project_directory, "windows/${index_array}.vcf", sep=""), file=a.script, append=T)
	write("", file=a.script, append=T)
	write("done", file=a.script, append=T)
