

diplo_vcf_list <- scan("diplo_vcf_list.txt", what="character")

diplo_vcf_list1 <- diplo_vcf_list[seq(from=1, to=(length(diplo_vcf_list) / 2), by=1) * 2 - 1]
diplo_vcf_list2 <- diplo_vcf_list[seq(from=1, to=(length(diplo_vcf_list) / 2), by=1) * 2]

# loop for each chromosome in the analysis
for(a in 1:length(diplo_vcf_list1)) {
	
	# albescens
	a_vcf <- paste("vcf/", diplo_vcf_list2[a], sep="")
	a_diplo_name <- sapply(strsplit(diplo_vcf_list2[a], "\\."), "[[", 1)
	a_script_name <- paste("a_", a_diplo_name, ".sh", sep="")
	a_chrom <- diplo_vcf_list1[a]
	a_output_name <- paste(a_chrom, "-alb.predictions_50", sep="")
	a_fvec <- paste("fvec_50kbp/", a_diplo_name, ".fvec", sep="")
	
	# write
	write("#!/bin/sh", file= a_script_name, ncolumns=1)
	write("#$ -V", file= a_script_name, ncolumns=1, append=T)
	write("#$ -cwd", file= a_script_name, ncolumns=1, append=T)
	write("#$ -S /bin/bash", file= a_script_name, ncolumns=1, append=T)
	write(paste("#$ -N aa_", a_diplo_name, sep=""), file= a_script_name, ncolumns=1, append=T)
	write("#$ -q omni", file= a_script_name, ncolumns=1, append=T)
	write("#$ -pe sm 4", file= a_script_name, ncolumns=1, append=T)
	write("#$ -P quanah", file= a_script_name, ncolumns=1, append=T)
	write("#$ -l h_rt=48:00:00", file= a_script_name, ncolumns=1, append=T)
	write("#$ -l h_vmem=8G", file= a_script_name, ncolumns=1, append=T)
	write("", file= a_script_name, ncolumns=1, append=T)
	write("source activate selection2", file= a_script_name, ncolumns=1, append=T)
	write("", file= a_script_name, ncolumns=1, append=T)
	write("python ~/diploSHIC/diploSHIC.py predict \\", file= a_script_name, ncolumns=1, append=T)
	write("albescensModel.json albescensModel.weights.hdf5 \\", file= a_script_name, ncolumns=1, append=T)
	write(paste(a_fvec, " ",  a_output_name, sep=""), file= a_script_name, ncolumns=1, append=T)
	write("", file= a_script_name, ncolumns=1, append=T)

	a_output_name <- paste(a_chrom, "-alb.predictions_20", sep="")
	a_fvec <- paste("fvec_20kbp/", a_diplo_name, ".fvec", sep="")
	write("", file= a_script_name, ncolumns=1, append=T)
	write("", file= a_script_name, ncolumns=1, append=T)
	write("python ~/diploSHIC/diploSHIC.py predict \\", file= a_script_name, ncolumns=1, append=T)
	write("albescensModel.json albescensModel.weights.hdf5 \\", file= a_script_name, ncolumns=1, append=T)
	write(paste(a_fvec, " ",  a_output_name, sep=""), file= a_script_name, ncolumns=1, append=T)
	write("", file= a_script_name, ncolumns=1, append=T)
		
	
	# americana
	a_vcf <- paste("vcf/", diplo_vcf_list2[a], sep="")
	a_vcf <- gsub("albescens", "americana", a_vcf)
	a_diplo_name <- sapply(strsplit(diplo_vcf_list2[a], "\\."), "[[", 1)
	a_diplo_name <- gsub("albescens", "americana", a_diplo_name)
	a_script_name <- paste("a_", a_diplo_name, ".sh", sep="")
	a_chrom <- diplo_vcf_list1[a]
	a_output_name <- paste(a_chrom, "-ame.predictions_50", sep="")
	a_fvec <- paste("fvec_50kbp/", a_diplo_name, ".fvec", sep="")
	
	# write
	write("#!/bin/sh", file= a_script_name, ncolumns=1)
	write("#$ -V", file= a_script_name, ncolumns=1, append=T)
	write("#$ -cwd", file= a_script_name, ncolumns=1, append=T)
	write("#$ -S /bin/bash", file= a_script_name, ncolumns=1, append=T)
	write(paste("#$ -N aa_", a_diplo_name, sep=""), file= a_script_name, ncolumns=1, append=T)
	write("#$ -q omni", file= a_script_name, ncolumns=1, append=T)
	write("#$ -pe sm 4", file= a_script_name, ncolumns=1, append=T)
	write("#$ -P quanah", file= a_script_name, ncolumns=1, append=T)
	write("#$ -l h_rt=48:00:00", file= a_script_name, ncolumns=1, append=T)
	write("#$ -l h_vmem=8G", file= a_script_name, ncolumns=1, append=T)
	write("", file= a_script_name, ncolumns=1, append=T)
	write("source activate selection2", file= a_script_name, ncolumns=1, append=T)
	write("", file= a_script_name, ncolumns=1, append=T)
	write("python ~/diploSHIC/diploSHIC.py predict \\", file= a_script_name, ncolumns=1, append=T)
	write("americanaModel.json americanaModel.weights.hdf5 \\", file= a_script_name, ncolumns=1, append=T)
	write(paste(a_fvec, " ",  a_output_name, sep=""), file= a_script_name, ncolumns=1, append=T)
	write("", file= a_script_name, ncolumns=1, append=T)

	a_output_name <- paste(a_chrom, "-ame.predictions_20", sep="")
	a_fvec <- paste("fvec_20kbp/", a_diplo_name, ".fvec", sep="")
	write("", file= a_script_name, ncolumns=1, append=T)
	write("", file= a_script_name, ncolumns=1, append=T)
	write("python ~/diploSHIC/diploSHIC.py predict \\", file= a_script_name, ncolumns=1, append=T)
	write("americanaModel.json americanaModel.weights.hdf5 \\", file= a_script_name, ncolumns=1, append=T)
	write(paste(a_fvec, " ",  a_output_name, sep=""), file= a_script_name, ncolumns=1, append=T)
	write("", file= a_script_name, ncolumns=1, append=T)
	
}


