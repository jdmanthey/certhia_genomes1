# concatenates the prediction output files

x20_alb <- list.files("predictions", pattern=("*alb.predictions_20"), full.names=T)
x20_ame <- list.files("predictions", pattern=("*ame.predictions_20"), full.names=T)
x50_alb <- list.files("predictions", pattern=("*alb.predictions_50"), full.names=T)
x50_ame <- list.files("predictions", pattern=("*ame.predictions_50"), full.names=T)

x20_alb_output <- c()
x20_ame_output <- c()
x50_alb_output <- c()
x50_ame_output <- c()
for(a in 1:length(x20_alb)) {
	a_rep <- read.table(x20_alb[a], sep="\t", header=T)
	x20_alb_output <- rbind(x20_alb_output, a_rep)
	a_rep <- read.table(x20_ame[a], sep="\t", header=T)
	x20_ame_output <- rbind(x20_ame_output, a_rep)
	a_rep <- read.table(x50_alb[a], sep="\t", header=T)
	x50_alb_output <- rbind(x50_alb_output, a_rep)
	a_rep <- read.table(x50_ame[a], sep="\t", header=T)
	x50_ame_output <- rbind(x50_ame_output, a_rep)
}


write.table(x20_alb_output, file="alb_predictions_20kbp.txt", sep="\t", quote=F, row.names=F)
write.table(x20_ame_output, file="ame_predictions_20kbp.txt", sep="\t", quote=F, row.names=F)
write.table(x50_alb_output, file="alb_predictions_50kbp.txt", sep="\t", quote=F, row.names=F)
write.table(x50_ame_output, file="ame_predictions_50kbp.txt", sep="\t", quote=F, row.names=F)
