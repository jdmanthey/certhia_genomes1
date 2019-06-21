x_files <- list.files(pattern="*msOut")
x_names <- substr(x_files, 1, nchar(x_files) - 6)
x_names <- gsub("alb", "albescens", x_names)
x_names <- gsub("ame", "americana", x_names)
x_names <- paste(x_names, ".txt", sep="")

num_individuals_albescens <- 14
num_individuals_americana <- 24

# we will keep the first 1000 simulations that finished for each 
# loop for each file
for(a in 1:length(x_files)) {
	# read in file
	con <- file(x_files[a], open="r")
	a_rep <- readLines(con)
	
	# determine which species this file is for
	if(length(grep("albescens", x_names[a])) == 1) {
		a_number <- num_individuals_albescens + 1
	} else {
		a_number <- num_individuals_americana + 1		
	}
	
	# find the lines of the beginning of each simulated sequence block
	a_lines <- seq(from=1, to=length(a_rep), by=1)[grepl("segsites", a_rep)]
	
	# remove simulations with zero segregating sites from the list
	a_lines <- a_lines[a_rep[a_lines] != "segsites: 0"]
	print(length(a_lines))
	
	# write header
	write(a_rep[1], file=x_names[a], ncolumns=1)
	write(a_rep[2], file=x_names[a], ncolumns=1, append=T)
	
	# get 1000 simulations for each file
	counter <- 1
	while(counter <= 1000) {
		write("", file=x_names[a], ncolumns=1, append=T)
		write("//", file=x_names[a], ncolumns=1, append=T)
		
		# get this simulation info
		b_rep <- a_rep[seq(from=a_lines[1], to=(a_lines[1] + a_number), by=1)]
		# write each of the lines
		for(b in 1:length(b_rep)) {
			write(b_rep[b], file=x_names[a], ncolumns=1, append=T)
		}
		
		# remove that record from a_lines object
		a_lines <- a_lines[2:length(a_lines)]
		
		# add to counter 
		counter <- counter + 1
		
	}
	
	
}
