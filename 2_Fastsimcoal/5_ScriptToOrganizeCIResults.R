#Script to analyse fastsimcoal results (after parametric Bootstrap - Andrea Thomaz, Feb 28, 2017)

require(plyr)
##Reading the bestlhoods files and creating a dataframe
wd <- "~/Dropbox/Comparative/Fastsimcoal/input_boot_otherModels/Mima_SOUTH_IM/Mima_SOUTH_IM_boot/"
setwd(wd)
best <- list.files(wd, pattern = ".bestlhoods")
best_DF <- data.frame(ANCSIZE = numeric(),NPOP=numeric(),TDIV=numeric(),
                      MaxEstLhood=numeric(),MaxObsLhood=numeric())
input_name <- "Mima_SOUTH_IM_"
simple_name <- gsub(paste(input_name,"|.bestlhoods", sep=""), "", best)
name_split <- strsplit(simple_name, "_")
name_df <- data.frame(t(matrix(unlist(name_split), nrow=2)))
name_df <- cbind(gsub(".bestlhoods", "", best), name_df)
colnames(name_df) <- c("file_name","MSFS", "rep")

for (i in best){
  temp <- read.table(i, header = T)
  best_DF <- rbind(best_DF, temp)
}
#Table with all simulations
best_DF <- cbind(name_df, best_DF)
#write.table(best_DF,paste("../",input_name,"Boot.txt", sep=""),sep="\t", row.names = FALSE)

quantile(best_DF$ANCSIZE,c(.025,.975))
quantile(best_DF$NPOP,c(.025,.975))


#test_max <- subset(best_DF, MSFS == 100)
#table with the best likelihood replicate of each simuation
best_like <- ddply(best_DF, ~MSFS, function(x){
  x[which.max(x$MaxEstLhood),]})
#write.table(best_like,paste("../",input_name,"FinalBest.txt", sep=""),sep="\t", row.names = FALSE)

quantile(best_like$ANCSIZE,c(.025,.975))
hist(best_like$ANCSIZE)
quantile(best_like$NPOP,c(.025,.975))
hist(best_like$NPOP)
quantile(best_like$TDIV,c(.025,.975))
hist(best_like$TDIV)
quantile(best_like$MIG12,c(.025,.975))
quantile(best_like$MIG21,c(.025,.975))


#Creating table to be used with script Cal95CI.py (and calculatw 95% CI)
##SC
Table_for_CI<- best_like[,c("MSFS","ANCSIZE","NPOP","TDIV","TAdm","A12","A21","MaxEstLhood","MaxObsLhood")]
##IM
Table_for_CI<- best_like[,c("MSFS","ANCSIZE","NPOP","TDIV","MIG12","MIG21","MaxEstLhood","MaxObsLhood")]

head(Table_for_CI)
write.table(Table_for_CI,paste("../",input_name,"bestlhoods.summary", sep=""),sep="\t", row.names = FALSE)
