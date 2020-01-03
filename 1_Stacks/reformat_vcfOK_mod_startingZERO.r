#November 17, 2015 - Andrea Thomaz
#read .vcf file from stacks1.41 output (WITHOUT write_random_SNP flag) for:
#plot the frequency of variable sites per position along all loci
#calculate theta based on number of segregating sites and individuals to create blacklist to delete very variable loci.
#note that this script will not work for stacks2 or above.

require(plyr)
require(pegas)

#READ VCF
setwd('/Users/andreathomaz/Desktop/Comparative/H_boulengeri/stacks/')
data <- read.table('batch_1.vcf', header = FALSE, sep = "\t")
head(data,10)

#SEQUENCE LENGTH
seq_len <- 140 #MODIFY HERE according the length of the sequence

#creates dataframe with loci ID, the variable positions and the number of individuals in each loci
new_data <- data.frame(loci_ID = data[,3], pos_vcf = data[,2],
                       pos = (data[,2] - seq_len*(data[,3]-1))-2, 
                       ind = rowSums(data[,10:length(data)] != "./.:0:.,.:.,.,."))
head(new_data)
min(new_data$pos)
max(new_data$pos)
length(unique(new_data$loci_ID))

#saving graph with frequency of variable sites along the loci
#pdf("./SNPdistr_pos140bp.pdf")
hist(new_data$pos, xlim = c(-1,seq_len), breaks = c(seq(-1, seq_len-1, by=1)), xlab = 'Position along the loci', main = 'The position of segregating sites');
abline(2050, 0, col = "red")
abline(v = 120, col = "red")
#dev.off()

#BASE ON THE GRAPH, CHOOSE HOW MANY POSITION TO DELETE FROM THE END
to_del <- 20 #how many sites to delete in the end of the sequence
seq_len_cut <- seq_len - to_del

whitelist <- subset(new_data, pos < seq_len_cut)[,c(1,3,4)]
white_temp <- whitelist$pos - 4
#pdf("./SNPdistr_pos_cutto120bp.pdf")
hist(white_temp, xlim = c(1,seq_len_cut), breaks = c(seq(-1, seq_len_cut -5 , by=1)), xlab = 'Position along the loci', main = 'The position of segregating sites');
#dev.off()

#calculating theta
var.sites <- count(whitelist, "loci_ID")
length(var.sites$loci_ID)
max(var.sites$freq)
theta_calc <- merge(unique(whitelist[,-2]), var.sites, by = "loci_ID")
theta_calc$theta <- 0
for (i in 1:length(theta_calc$theta)){
  theta_calc[i,4] <- theta.s(theta_calc$freq[i], theta_calc$ind[i])/seq_len_cut
}

#calculating the 95% quantile to exclude loci extremely variable
quant <- quantile(theta_calc$theta, probs=0.95) #set the value to be 
quant
#pdf("./theta120bp.pdf")
hist(theta_calc$theta)
abline(v = quant, col="red")
#dev.off()

#what is the maximum number of mutations in a loci
x <- subset(theta_calc, theta < quant)
max(x$freq)

#saving whitelist for re-run populations in stacks
blacklist <- subset(theta_calc, theta > quant)[,1]
#write.table(blacklist, file="blacklist.txt", sep = '\n', row.names = F, col.names = F)

whitelist$blacklist <- match(whitelist$loci_ID, blacklist, nomatch = 0)
whitelist_final <- subset(whitelist, blacklist == 0)
length(unique(whitelist_final$loci_ID))
length(whitelist_final$loci_ID)
write.table(whitelist_final[,1:2], file="whitelist_Hboulengeri.txt", sep = '\t', row.names = F, col.names = F)

