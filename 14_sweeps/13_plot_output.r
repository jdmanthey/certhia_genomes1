options(scipen=999)
require(stats)

x1 <- read.table("alb_predictions_50kbp.txt", sep="\t", header=T, stringsAsFactors=F)
x2 <- read.table("ame_predictions_50kbp.txt", sep="\t", header=T, stringsAsFactors=F)

# window size:
kbp_windows <- 50000

# chroms 
chroms <- unique(x1[,1])
# reorganize the chromosome order
chroms <- chroms[c(2,5,1,3,6,14,7,9,8,10,11,15,12,13,16,18,19,21,20,23,26,17,24,22,4)]

# reorganize chromosome order in matrices
new_x1 <- c()
new_x2 <- c()
for(a in 1:length(chroms)) {
	new_x1 <- rbind(new_x1, x1[x1[,1] == chroms[a], ])
	new_x2 <- rbind(new_x2, x2[x2[,1] == chroms[a], ])
}
# put back into main matrix names
x1 <- new_x1
x2 <- new_x2
# rename rownames
rownames(x1) <- seq(from=1, to=nrow(x1), by=1)
rownames(x2) <- seq(from=1, to=nrow(x2), by=1)

# add to counts of the window locations for plotting
count <- 0
for(a in 1:length(chroms)) {
	# length of this chromosome 
	a_rep_max <- max(x2[x2[,1] == chroms[a],3])
	
	# add the distance already covered by count
	x1[x1[,1] == chroms[a], 2] <- x1[x1[,1] == chroms[a], 2] + count
	x1[x1[,1] == chroms[a], 3] <- x1[x1[,1] == chroms[a], 3] + count
	x2[x2[,1] == chroms[a], 2] <- x2[x2[,1] == chroms[a], 2] + count
	x2[x2[,1] == chroms[a], 3] <- x2[x2[,1] == chroms[a], 3] + count
	
	# add length of this chromosome to the count
	count <- count + a_rep_max
}

#check that count addition worked
tail(x1)
tail(x2)



# plot pie charts
par(mfrow=c(1,2))
par(mar=c(1,1,1,1))
pie(table(x2$predClass), labels=NA)
pie(table(x1$predClass), labels=NA)




# where are the shared hard sweeps 

hard1 <- x1[x1$predClass == "hard" | x1$predClass == "linkedHard", ]
locations_hard1 <- paste(hard1[,1], hard1[,2], sep="_")
hard2 <- x2[x2$predClass == "hard" | x2$predClass == "linkedHard", ]
locations_hard2 <- paste(hard2[,1], hard2[,2], sep="_")
# shared hard sweep or linked hard
hard1[locations_hard1 %in% locations_hard2,]

# how many windows are neutral in both lineages?
neutral1 <- x1[x1$predClass == "neutral", ]
neutral2 <- x2[x2$predClass == "neutral", ]
locations_neutral1 <- paste(neutral1[,1], neutral1[,2], sep="_")
locations_neutral2 <- paste(neutral2[,1], neutral2[,2], sep="_")
table(locations_neutral1 %in% locations_neutral2)




# plot neutrality and classifications


five_class_colors <- c("lightgray", "lightblue", "lightcoral", "navyblue", "red3")


window_size <- 10


# what are the unique chromosomes and their bounding areas for plotting?
chr <- unique(x1[,1])
chr_polygons_1 <- list()
chr_polygons_2 <- list()
# make the plotting polygons
for(a in 1:length(chr)) {
	a1 <- rownames(x1)[x1[,1] == chr[a]]
	a2 <- a1[length(a1)]
	a1 <- a1[1]
	chr_polygons_1[[a]] <- rbind(c(a1, 0), c(a2, 0), c(a2, 1), c(a1, 1), c(a1, 0))
	a1 <- rownames(x2)[x2[,1] == chr[a]]
	a2 <- a1[length(a1)]
	a1 <- a1[1]
	chr_polygons_2[[a]] <- rbind(c(a1, 0), c(a2, 0), c(a2, 1), c(a1, 1), c(a1, 0))
}





# set up plotting dimensions
par(mfrow=c(4,1))
par(mar=c(0.5,5,1,0))


# plot neutrality americana
plot(c(-1,-1), ylim=c(0,1), xlim=c(1, nrow(x2)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="Prob. Neutral")
odd <- 0
for(a in 1:length(chr_polygons_2)) {
	if(odd == 1) {
		polygon(chr_polygons_2[[a]], col="snow2", border="white")
		odd <- 0	
	} else {
		odd <- 1
	}
}
# plot
points(rownames(x2), x2$prob.neutral., pch=19, cex=0.1, col="sienna1")	

# sliding windows
total_rep <- c()
place_rep <- c()
sliding_windows <- ceiling(as.numeric(rownames(x2)[nrow(x2)]) / window_size)
location <- window_size/2
for(b in 1:sliding_windows) {
	
	# get center location
	b_rep <- x2[location,2]
	# take surrounding lines in matrix
	lines_to_take <- (location-(window_size/2 - 1)):(location+(window_size/2))
	lines_to_take <- lines_to_take[lines_to_take >= 1 & lines_to_take <= nrow(x2)]
	b_rep2 <- x2[lines_to_take,]
	# see if they are within the needed window size
	b_rep2 <- b_rep2[b_rep2[,2] >= (b_rep - kbp_windows * (window_size/2)) & b_rep2[,2] <= (b_rep + kbp_windows * (window_size/2)),]
	
	# get values of window
	total_rep <- c(total_rep, mean(b_rep2$prob.neutral.))
	# get value for x axis
	place_rep <- c(place_rep, location)
	# add to row locator
	location <- location + window_size

}
lines(place_rep, total_rep, lwd=0.6, col="sienna4")

# plot classification of window
plot(c(-1,-1), ylim=c(0,1), xlim=c(1, max(x2[,3])), xaxt="n", yaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="Window Classification")

for(a in 1:nrow(x2)) {
	if(x2$predClass[a] == "neutral") {
		plot_color <- five_class_colors[1]
		line_points <- rbind(c(x2$classifiedWinStart[a], 0), c(x2$classifiedWinStart[a], 0.2))
	} else if(x2$predClass[a] == "linkedSoft") {
		plot_color <- five_class_colors[2]
		line_points <- rbind(c(x2$classifiedWinStart[a], 0.2), c(x2$classifiedWinStart[a], 0.4))
	} else if(x2$predClass[a] == "linkedHard") {
		plot_color <- five_class_colors[3]
		line_points <- rbind(c(x2$classifiedWinStart[a], 0.6), c(x2$classifiedWinStart[a], 0.8))
	} else if(x2$predClass[a] == "soft") {
		plot_color <- five_class_colors[4]
		line_points <- rbind(c(x2$classifiedWinStart[a], 0.4), c(x2$classifiedWinStart[a], 0.6))
	} else {
		plot_color <- five_class_colors[5]
		line_points <- rbind(c(x2$classifiedWinStart[a], 0.8), c(x2$classifiedWinStart[a], 1))
	}
	
	lines(line_points, col=plot_color, lwd=0.08)
}



# plot neutrality albescens
plot(c(-1,-1), ylim=c(0,1), xlim=c(1, nrow(x1)), xaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="Prob. Neutral")
odd <- 0
for(a in 1:length(chr_polygons_2)) {
	if(odd == 1) {
		polygon(chr_polygons_2[[a]], col="snow2", border="white")
		odd <- 0	
	} else {
		odd <- 1
	}
}
# plot
points(rownames(x1), x1$prob.neutral., pch=19, cex=0.1, col="sienna1")	

# sliding windows
total_rep <- c()
place_rep <- c()
sliding_windows <- ceiling(as.numeric(rownames(x1)[nrow(x1)]) / window_size)
location <- window_size/2
for(b in 1:sliding_windows) {
	
	# get center location
	b_rep <- x1[location,2]
	# take surrounding lines in matrix
	lines_to_take <- (location-(window_size/2 - 1)):(location+(window_size/2))
	lines_to_take <- lines_to_take[lines_to_take >= 1 & lines_to_take <= nrow(x1)]
	b_rep2 <- x1[lines_to_take,]
	# see if they are within the needed window size
	b_rep2 <- b_rep2[b_rep2[,2] >= (b_rep - kbp_windows * (window_size/2)) & b_rep2[,2] <= (b_rep + kbp_windows * (window_size/2)),]
	
	# get values of window
	total_rep <- c(total_rep, mean(b_rep2$prob.neutral.))
	# get value for x axis
	place_rep <- c(place_rep, location)
	# add to row locator
	location <- location + window_size

}
lines(place_rep, total_rep, lwd=0.6, col="sienna4")

# plot classification of window
plot(c(-1,-1), ylim=c(0,1), xlim=c(1, max(x1[,3])), xaxt="n", yaxt="n", col="white", bty="n", cex.axis=1.1, cex.lab=1.3, ylab="Window Classification")

for(a in 1:nrow(x1)) {
	if(x1$predClass[a] == "neutral") {
		plot_color <- five_class_colors[1]
		line_points <- rbind(c(x1$classifiedWinStart[a], 0), c(x1$classifiedWinStart[a], 0.2))
	} else if(x1$predClass[a] == "linkedSoft") {
		plot_color <- five_class_colors[2]
		line_points <- rbind(c(x1$classifiedWinStart[a], 0.2), c(x1$classifiedWinStart[a], 0.4))
	} else if(x1$predClass[a] == "linkedHard") {
		plot_color <- five_class_colors[3]
		line_points <- rbind(c(x1$classifiedWinStart[a], 0.6), c(x1$classifiedWinStart[a], 0.8))
	} else if(x1$predClass[a] == "soft") {
		plot_color <- five_class_colors[4]
		line_points <- rbind(c(x1$classifiedWinStart[a], 0.4), c(x1$classifiedWinStart[a], 0.6))
	} else {
		plot_color <- five_class_colors[5]
		line_points <- rbind(c(x1$classifiedWinStart[a], 0.8), c(x1$classifiedWinStart[a], 1))
	}
	
	lines(line_points, col=plot_color, lwd=0.08)
}


