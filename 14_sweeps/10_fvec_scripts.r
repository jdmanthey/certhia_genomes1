

diplo_vcf_list <- scan("diplo_vcf_list.txt", what="character")
fasta_index <- read.table("06_certhia_reordered.fasta.fai", sep="\t", stringsAsFactors=F)

diplo_vcf_list1 <- diplo_vcf_list[seq(from=1, to=(length(diplo_vcf_list) / 2), by=1) * 2 - 1]
diplo_vcf_list2 <- diplo_vcf_list[seq(from=1, to=(length(diplo_vcf_list) / 2), by=1) * 2]

# loop for each chromosome in the analysis
for(a in 1:length(diplo_vcf_list1)) {
	
	# albescens
	a_vcf <- paste("vcf/", diplo_vcf_list2[a], sep="")
	a_diplo_name <- sapply(strsplit(diplo_vcf_list2[a], "\\."), "[[", 1)
	a_output_name <- paste("a_", a_diplo_name, ".sh", sep="")
	a_chrom <- diplo_vcf_list1[a]
	a_length <- fasta_index[fasta_index[,1] == a_chrom,2]
	a_fvec <- paste("fvec_50kbp/", a_diplo_name, ".fvec", sep="")
	
	# write
	write("#!/bin/sh", file=a_output_name, ncolumns=1)
	write("#$ -V", file=a_output_name, ncolumns=1, append=T)
	write("#$ -cwd", file=a_output_name, ncolumns=1, append=T)
	write("#$ -S /bin/bash", file=a_output_name, ncolumns=1, append=T)
	write(paste("#$ -N a_", a_diplo_name, sep=""), file=a_output_name, ncolumns=1, append=T)
	write("#$ -q omni", file=a_output_name, ncolumns=1, append=T)
	write("#$ -pe sm 4", file=a_output_name, ncolumns=1, append=T)
	write("#$ -P quanah", file=a_output_name, ncolumns=1, append=T)
	write("#$ -l h_rt=48:00:00", file=a_output_name, ncolumns=1, append=T)
	write("#$ -l h_vmem=8G", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)
	write("source activate selection2", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)
	write("python ~/diploSHIC/diploSHIC.py fvecVcf diploid \\", file=a_output_name, ncolumns=1, append=T)
	write(paste(a_vcf, " ", a_chrom, " ", a_length, " \\", sep=""), file=a_output_name, ncolumns=1, append=T)
	write(paste(a_fvec, " --targetPop alb --sampleToPopFileName popmap.txt --winSize 550000 \\", sep=""), 
		file=a_output_name, ncolumns=1, append=T)
	write("--maskFileName diplo_mask.fasta.gz", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)

	# albescens second size
	a_fvec <- paste("fvec_20kbp/", a_diplo_name, ".fvec", sep="")
	
	# write
	write("", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)
	write("python ~/diploSHIC/diploSHIC.py fvecVcf diploid \\", file=a_output_name, ncolumns=1, append=T)
	write(paste(a_vcf, " ", a_chrom, " ", a_length, " \\", sep=""), file=a_output_name, ncolumns=1, append=T)
	write(paste(a_fvec, " --targetPop alb --sampleToPopFileName popmap.txt --winSize 220000 \\", sep=""), 
		file=a_output_name, ncolumns=1, append=T)
	write("--maskFileName diplo_mask.fasta.gz", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)

	
	
	
	# americana
	a_vcf <- paste("vcf/", diplo_vcf_list2[a], sep="")
	a_vcf <- gsub("albescens", "americana", a_vcf)
	a_diplo_name <- sapply(strsplit(diplo_vcf_list2[a], "\\."), "[[", 1)
	a_diplo_name <- gsub("albescens", "americana", a_diplo_name)
	a_output_name <- paste("a_", a_diplo_name, ".sh", sep="")
	a_chrom <- diplo_vcf_list1[a]
	a_length <- fasta_index[fasta_index[,1] == a_chrom,2]
	a_fvec <- paste("fvec_50kbp/", a_diplo_name, ".fvec", sep="")
	
	# write
	write("#!/bin/sh", file=a_output_name, ncolumns=1)
	write("#$ -V", file=a_output_name, ncolumns=1, append=T)
	write("#$ -cwd", file=a_output_name, ncolumns=1, append=T)
	write("#$ -S /bin/bash", file=a_output_name, ncolumns=1, append=T)
	write(paste("#$ -N a_", a_diplo_name, sep=""), file=a_output_name, ncolumns=1, append=T)
	write("#$ -q omni", file=a_output_name, ncolumns=1, append=T)
	write("#$ -pe sm 4", file=a_output_name, ncolumns=1, append=T)
	write("#$ -P quanah", file=a_output_name, ncolumns=1, append=T)
	write("#$ -l h_rt=48:00:00", file=a_output_name, ncolumns=1, append=T)
	write("#$ -l h_vmem=8G", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)
	write("source activate selection2", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)
	write("python ~/diploSHIC/diploSHIC.py fvecVcf diploid \\", file=a_output_name, ncolumns=1, append=T)
	write(paste(a_vcf, " ", a_chrom, " ", a_length, " \\", sep=""), file=a_output_name, ncolumns=1, append=T)
	write(paste(a_fvec, " --targetPop ame --sampleToPopFileName popmap.txt --winSize 550000 \\", sep=""), 
		file=a_output_name, ncolumns=1, append=T)
	write("--maskFileName diplo_mask.fasta.gz", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)

	# americana second size
	a_fvec <- paste("fvec_20kbp/", a_diplo_name, ".fvec", sep="")
	
	# write
	write("", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)
	write("python ~/diploSHIC/diploSHIC.py fvecVcf diploid \\", file=a_output_name, ncolumns=1, append=T)
	write(paste(a_vcf, " ", a_chrom, " ", a_length, " \\", sep=""), file=a_output_name, ncolumns=1, append=T)
	write(paste(a_fvec, " --targetPop ame --sampleToPopFileName popmap.txt --winSize 220000 \\", sep=""), 
		file=a_output_name, ncolumns=1, append=T)
	write("--maskFileName diplo_mask.fasta.gz", file=a_output_name, ncolumns=1, append=T)
	write("", file=a_output_name, ncolumns=1, append=T)

}





