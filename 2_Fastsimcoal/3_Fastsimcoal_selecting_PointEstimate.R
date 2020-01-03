#Script to analyse fastsimcoal results - Andrea Thomaz, Mar 2017
library(plyr)

##Needs to have all .bestlhoods in the same folder and 
wd <- "~/Dropbox/Comparative/Fastsimcoal/input_fastsimcoal/SecContact_model/Results/"
setwd(wd)
filez <- list.files(wd)
filez

best_DF <- data.frame(SP = character(), FILE = character(), ANCSIZE = numeric(),NPOP=numeric(),TDIV=numeric(),
                      MaxEstLhood=numeric(),MaxObsLhood=numeric())

for (i in filez){
    temp <- read.table(i, header = T)

    temp <- cbind(SP = gsub("[_][0-9]+","",gsub(".bestlhoods", "", i)),
                  FILE = gsub("[A-Za-z_]","",gsub(".bestlhoods", "", i)), temp)
    best_DF <- rbind(best_DF, temp)
  }

write.csv(best_DF, "../Fastsimcoal_AllRunsPerModel.csv")


res <- ddply(best_DF, ~SP, function(x){
             data.frame(ANCSIZE_mean = mean(x$ANCSIZE), ANCSIZE_sd = sd(x$ANCSIZE),
             NPOP_mean = mean(x$NPOP), NPOP_sd = sd(x$NPOP),
             TDIV_mean = mean(x$TDIV), TDIV_sd = sd(x$TDIV))
             })

res <- ddply(best_DF, ~SP, function(x){
              x[which.max(x$MaxEstLhood),]})
write.csv(res, "../Fatsimcoal_PointEstimate.csv")
