output_directory <- "cds_fasta_files"
output_summary <- "cds_summary_and_mapping.txt"
dir.create(output_directory)

# read in fasta files
certhia_fasta <- scan("certhia_cds_renamed.fasta", what="character")
certhia_fasta <- gsub(">", "", certhia_fasta)
ficedula_fasta <- scan("ficedula_cds_renamed.fasta", what="character")
ficedula_fasta <- gsub(">", "", ficedula_fasta)
parus_fasta <- scan("parus_cds_renamed.fasta", what="character")
parus_fasta <- gsub(">", "", parus_fasta)
taeniopygia_fasta <- scan("taeniopygia_cds_renamed.fasta", what="character")
taeniopygia_fasta <- gsub(">", "", taeniopygia_fasta)

# put each fasta into a dataframe with column 1 = name and column 2 = sequence
certhia_cds <- data.frame(name=as.character(certhia_fasta[1:(length(certhia_fasta) / 2) * 2 - 1]),
							sequence=as.character(certhia_fasta[1:(length(certhia_fasta) / 2) * 2]))
ficedula_cds <- data.frame(name=as.character(ficedula_fasta[1:(length(ficedula_fasta) / 2) * 2 - 1]),
							sequence=as.character(ficedula_fasta[1:(length(ficedula_fasta) / 2) * 2]))
parus_cds <- data.frame(name=as.character(parus_fasta[1:(length(parus_fasta) / 2) * 2 - 1]),
							sequence=as.character(parus_fasta[1:(length(parus_fasta) / 2) * 2]))
taeniopygia_cds <- data.frame(name=as.character(taeniopygia_fasta[1:(length(taeniopygia_fasta) / 2) * 2 - 1]),
							sequence=as.character(taeniopygia_fasta[1:(length(taeniopygia_fasta) / 2) * 2]))
							


# read in blast files
q_certhia_s_ficedula <- read.table("q_certhia_s_ficedula.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_certhia_s_parus <- read.table("q_certhia_s_parus.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_certhia_s_taeniopygia <- read.table("q_certhia_s_taeniopygia.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_ficedula_s_certhia <- read.table("q_ficedula_s_certhia.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_parus_s_certhia <- read.table("q_parus_s_certhia.blast", sep="\t", stringsAsFactors=F, quote="\"")
q_taeniopygia_s_certhia <- read.table("q_taeniopygia_s_certhia.blast", sep="\t", stringsAsFactors=F, quote="\"")


# write header of summary file
write(c("certhia_gene_number", "ficedula_gene_name", "parus_gene_name", "taeniopygia_gene_name", "number"), file=output_summary, ncolumns=5)

# list each certhia predicted gene
certhia_genes <- as.character(unique(certhia_cds[,1]))

# start counter
counter <- 1

# loop for each predicted gene
for(a in 1:length(certhia_genes)) {
	
	# matches where the subject is the other species
	s_ficedula <- q_certhia_s_ficedula[grep(certhia_genes[a], q_certhia_s_ficedula[,1]),]
	s_parus <- q_certhia_s_parus[grep(certhia_genes[a], q_certhia_s_parus[,1]),]
	s_taeniopygia <- q_certhia_s_taeniopygia[grep(certhia_genes[a], q_certhia_s_taeniopygia[,1]),]
	# keep going if there is a match for all three species
	if(nrow(s_ficedula) > 0 & nrow(s_parus) > 0 & nrow(s_taeniopygia) > 0) {
		# check how many unique gene matches for each search
		if(length(unique(s_ficedula[,2])) == 1) {
			s_ficedula <- as.character(s_ficedula[1,2])
		} else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
			s_ficedula <- s_ficedula[s_ficedula[,11] == min(s_ficedula[,11]),]
			s_ficedula <- s_ficedula[s_ficedula[,12] == max(s_ficedula[,12]),]
			s_ficedula <- as.character(s_ficedula[1,2])
		}
		if(length(unique(s_parus[,2])) == 1) {
			s_parus <- as.character(s_parus[1,2])
		} else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
			s_parus <- s_parus[s_parus[,11] == min(s_parus[,11]),]
			s_parus <- s_parus[s_parus[,12] == max(s_parus[,12]),]
			s_parus <- as.character(s_parus[1,2])
		}
		if(length(unique(s_taeniopygia[,2])) == 1) {
			s_taeniopygia <- as.character(s_taeniopygia[1,2])
		} else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
			s_taeniopygia <- s_taeniopygia[s_taeniopygia[,11] == min(s_taeniopygia[,11]),]
			s_taeniopygia <- s_taeniopygia[s_taeniopygia[,12] == max(s_taeniopygia[,12]),]
			s_taeniopygia <- as.character(s_taeniopygia[1,2])
		}
		
		# search for the matched genes in the searches where certhia was the subject to check for reciprocal matches
		q_ficedula <- q_ficedula_s_certhia[grep(s_ficedula, q_ficedula_s_certhia[,1]),]
		q_parus <- q_parus_s_certhia[grep(s_parus, q_parus_s_certhia[,1]),]
		q_taeniopygia <- q_taeniopygia_s_certhia[grep(s_taeniopygia, q_taeniopygia_s_certhia[,1]),]
		# keep going if there is a match for all three species
		if(nrow(q_ficedula) > 0 & nrow(q_parus) > 0 & nrow(q_taeniopygia) > 0) {
			if(length(unique(q_ficedula[,2])) == 1) {
				q_ficedula <- as.character(q_ficedula[1,2])
			} else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
				q_ficedula <- q_ficedula[q_ficedula[,11] == min(q_ficedula[,11]),]
				q_ficedula <- q_ficedula[q_ficedula[,12] == max(q_ficedula[,12]),]
				q_ficedula <- as.character(q_ficedula[1,2])
			}
			if(length(unique(q_parus[,2])) == 1) {
				q_parus <- as.character(q_parus[1,2])
			} else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
				q_parus <- q_parus[q_parus[,11] == min(q_parus[,11]),]
				q_parus <- q_parus[q_parus[,12] == max(q_parus[,12]),]
				q_parus <- as.character(q_parus[1,2])
			}
			if(length(unique(q_taeniopygia[,2])) == 1) {
				q_taeniopygia <- as.character(q_taeniopygia[1,2])
			} else { # choose best match by minimum e-value, if two have same e-value, choose higher bit score
				q_taeniopygia <- q_taeniopygia[q_taeniopygia[,11] == min(q_taeniopygia[,11]),]
				q_taeniopygia <- q_taeniopygia[q_taeniopygia[,12] == max(q_taeniopygia[,12]),]
				q_taeniopygia <- as.character(q_taeniopygia[1,2])
			}
			
			# keep going if all three best matches 
			if(q_parus == certhia_genes[a] & q_taeniopygia == certhia_genes[a] & q_ficedula == certhia_genes[a]) {
				# write info to summary file
				write(c(certhia_genes[a], s_ficedula, s_parus, s_taeniopygia, counter), file=output_summary, ncolumns=5, append=T)
				
				# write the fasta files for each gene that got this far
				write(">certhia", file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1)
				write(as.character(certhia_cds[certhia_cds[,1] == certhia_genes[a],2]), 
					file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
				write(">ficedula", file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
				write(as.character(ficedula_cds[ficedula_cds[,1] == s_ficedula,2]), 
					file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
				write(">parus", file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
				write(as.character(parus_cds[parus_cds[,1] == s_parus,2]), 
					file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
				write(">taeniopygia", file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)
				write(as.character(taeniopygia_cds[taeniopygia_cds[,1] == s_taeniopygia,2]), 
					file=paste(output_directory, "/", counter, ".fasta", sep=""), ncolumns=1, append=T)	
									
				# add to counter
				counter <- counter + 1
			}
		}	
	}
}









