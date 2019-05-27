library(OmicCircos)
options(stringsAsFactors=FALSE)
library(pals)

x <- read.table("certhia_filtered.txt", sep="\t", header=T)
# rearrange the input file to the following order
# target name, start, end, query name, start, end, identity, length target, length query
x <- x[,c(10,1,2,11,3,4,7, 9,10)]

# subset scaffold names and start sites
x2 <- x[,c(1,2,4,5)]
# rename certhia scaffolds
x2[,3] <- paste("Ca_s", sapply(strsplit(x2[,3], "_"), "[[", 3), sep="")


# bin the segments of the circos plot
segf_lengths <- c()
x3 <- c()
for(a in 1:length(unique(x2[,1]))) {
	a_rep <- x2[x2[,1] == unique(x2[,1])[a], ]
	test <- unique(sort(a_rep[,2]))
	a_rep[,2] <- match(a_rep[,2], test)
	x3 <- rbind(x3, a_rep)
	
	segf_lengths <- rbind(segf_lengths, c(unique(x2[,1])[a], length(test)))
}


x4 <- c()
for(a in 1:length(unique(x2[,3]))) {
	a_rep <- x3[x3[,3] == unique(x2[,3])[a], ]
	test <- unique(sort(a_rep[,4]))
	a_rep[,4] <- match(a_rep[,4], test)
	x4 <- rbind(x4, a_rep)
	segf_lengths <- rbind(segf_lengths, c(unique(x2[,3])[a], length(test)))
}
segf_lengths <- data.frame(chrom=as.character(segf_lengths[,1]), bins=as.numeric(segf_lengths[,2]))


# create segments object
seg_f <- c()
for(a in 1:nrow(segf_lengths)) {
	a_name <- segf_lengths[a,1]
	a_lengths <- 0:(segf_lengths[a,2] - 1)
	a_lengths2 <- a_lengths + 1
	a_name <- rep(a_name, length(a_lengths))
	a_v <- rep("NA", length(a_lengths))
	a_note <- rep("NA", length(a_lengths))
	a_output <- cbind(a_name, a_lengths, a_lengths2, a_v, a_note)
	seg_f <- rbind(seg_f, a_output)
}
seg_f <- data.frame(seg.name=as.character(seg_f[,1]), seg.Start=as.numeric(seg_f[,2]), seg.End=as.numeric(seg_f[,3]), the.v=as.character(seg_f[,4]), Note=as.character(seg_f[,5]))


# make the links object in the right format
link_names <- paste("n", seq(from=1, to=nrow(x4), by=1), sep="")
link_v <- data.frame(seg1=as.character(x4[,1]), po1=as.numeric(x4[,2]), name1=as.character(link_names), seg2=as.character(x4[,3]), po2=as.numeric(x4[,4]), name2=as.character(link_names))
# subset for easier plotting
link_subset <- link_v[sample(1:nrow(link_v), floor(nrow(link_v) / 100)), ]


# make the angular database
seg_names <- sort(unique(segf_lengths[,1]))

# modify for particular plot if necessary:
seg_names <- c("Ca_s2", "Ca_s5", "Ca_s1", "Ca_s3", "Ca_s6", "Ca_s14", "Ca_s7", "Ca_s9", "Ca_s8", "Ca_s10", "Ca_s11", "Ca_s15", "Ca_s12", "Ca_s13", "Ca_s16", "Ca_s19", "Ca_s20", "Ca_s22", "Ca_s21", "Ca_s23", "Ca_s29", "Ca_s17", "Ca_s26", "Ca_s27", "Ca_s24", "Ca_s30", "Ca_s18", "Ca_s25", "Ca_s28", "Ca_s4", "Tg_Z", "Tg_28", "Tg_27", "Tg_26", "Tg_23", "Tg_22", "Tg_21", "Tg_20", "Tg_24", "Tg_19", "Tg_18", "Tg_17", "Tg_15", "Tg_14", "Tg_13", "Tg_12", "Tg_11", "Tg_10", "Tg_9", "Tg_8", "Tg_7", "Tg_6", "Tg_5", "Tg_4A", "Tg_4", "Tg_3", "Tg_2", "Tg_1A", "Tg_1") 
seg_colors <- c(rep("darkorange", 30), rep("navyblue", 29))

# maybe for plotting opposite?
#seg_names <- rev(seg_names)
#seg_colors <- rev(seg_colors)

db <- segAnglePo(seg_f, seg=seg_names)



# colors for plotting
cols <- polychrome(length(unique(link_subset[,4])))
cols_rgb <- col2rgb(cols)
new_cols <- c()
for(a in 1:ncol(cols_rgb)) {
	a_col <- rgb(cols_rgb[1,a], cols_rgb[2,a], cols_rgb[3,a], max=255, alpha=0.4*255)
	new_cols <- c(new_cols, a_col)
}
cols <- new_cols

plot_colors <- rep("un", nrow(link_subset))
for(a in 1:length(unique(link_subset[,4]))) {
	a_rep <- sort(unique(link_subset[,4]))[a]
	plot_colors[link_subset[,4] == a_rep] <- cols[a]
}



par(mar=c(2,2,2,2)) 
plot(c(1,800),  c(1,800), type="n", axes=FALSE, xlab="XLAB", ylab="YLAB", main="")
circos(R=360, cir=db, type="chr", col=seg_colors, print.chr.lab=TRUE, W=40, scale = F)
circos(R=350, cir=db, W=40, mapping=link_subset, type="link", lwd=0.75, col=plot_colors)























