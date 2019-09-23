
library(adegenet)



x <- read.structure("certhia_contact_10000b.stru",onerowperind=F,n.ind=24,n.loc=31948,ask=F,sep="\t")

grp <- find.clusters(x,max.n.clust=10,n.iter=1e5) #25 pc, 2 cluster
dapc1 <- dapc(x,grp$grp) # choose 10 pc 1 df
scatter(dapc1,scree.da=F,cex=1.4,clab=0)
compoplot(dapc1)





# read in structure results
x2 <- read.table("structure_results.txt", sep="\t", stringsAsFactors=F)

# sort by population
populations <- c("Utah", "Pinal", "Pinaleno", "SCatalina", "Chiricahua", "SRita", "Huachuca", "Mexico")

# make new matrix
output <- c()
for(a in 1:length(populations)) {
	a_rep1 <- x2[sapply(strsplit(x2[,1], "_"), "[[", 3) == populations[a], ]
	a_rep2 <- dapc1$posterior[sapply(strsplit(x2[,1], "_"), "[[", 3) == populations[a], ]
	output <- rbind(output, cbind(a_rep1, a_rep2))
}

par(mar=c(0.1, 0.1, 0.1, 0.1))
par(mfrow=c(2,1))
barplot(t(as.matrix(output[,2:3])), col=c("red", "blue"))
barplot(t(as.matrix(output[,4:5])), col=c("red", "blue"))
