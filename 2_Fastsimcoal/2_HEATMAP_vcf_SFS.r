library(ggplot2)
library(reshape2)
#library(RColorBrewer)

setwd("~/Dropbox/Comparative/Fastsimcoal/input_fastsimcoal/")

filez <- list.files()[c(1,3:8,10:12)]

#BrycoM
SFS <- readLines(paste("./",filez[1],"/",filez[1],"_MSFS.obs", sep=""))
sfs_num <- as.numeric(unlist(strsplit(SFS[2], "\t")))
sfs <- matrix(as.numeric(unlist(strsplit(SFS[3], "\t"))), nrow = sfs_num[2]+1, ncol = sfs_num[2]+1, byrow = T)
sfs10 <- log10(sfs)
sfs10[which(sfs10 == -Inf)] <- 0
rownames(sfs10) <- seq(0,sfs_num[2],1)
colnames(sfs10) <- seq(0,sfs_num[2],1)

BrycoM <- ggplot(melt(sfs10), aes(Var1,Var2, fill=value)) + 
                geom_raster() + 
                scale_fill_gradientn(colours = rev(terrain.colors(30)), limits = c(0,3.2))+
                scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
                xlab("") + ylab("")+
                labs(fill = "Alleles\n(log10)")+
                ggtitle("Bryconamericus - Central")+
                theme(legend.position = c(0.9, 0.65)) +
                guides(fill = guide_colourbar(barwidth = 1, barheight = 10)); BrycoM

#Bryco_SOUTH
SFS <- readLines(paste("./",filez[2],"/",filez[2],"_MSFS.obs", sep=""))
sfs_num <- as.numeric(unlist(strsplit(SFS[2], "\t")))
sfs <- matrix(as.numeric(unlist(strsplit(SFS[3], "\t"))), nrow = sfs_num[2]+1, ncol = sfs_num[2]+1, byrow = T)
sfs10 <- log10(sfs)
sfs10[which(sfs10 == -Inf)] <- 0
rownames(sfs10) <- seq(0,sfs_num[2],1)
colnames(sfs10) <- seq(0,sfs_num[2],1)

BrycoS <- ggplot(melt(sfs10), aes(Var1,Var2, fill=value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = rev(terrain.colors(30)), limits = c(0,3.2))+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  #xlab("Southern population") + ylab("Northern population")+
  xlab("") + ylab("")+
  ggtitle("Bryconamericus - South")+
  #labs(fill = "Alleles\n(log10)")+
  theme(legend.position = "none") +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 10)); BrycoS

#HboulM
SFS <- readLines(paste("./",filez[3],"/",filez[3],"_MSFS.obs", sep=""))
sfs_num <- as.numeric(unlist(strsplit(SFS[2], "\t")))
sfs <- matrix(as.numeric(unlist(strsplit(SFS[3], "\t"))), nrow = sfs_num[2]+1, ncol = sfs_num[2]+1, byrow = T)
sfs10 <- log10(sfs)
sfs10[which(sfs10 == -Inf)] <- 0
rownames(sfs10) <- seq(0,sfs_num[2],1)
colnames(sfs10) <- seq(0,sfs_num[2],1)

HboulM <- ggplot(melt(sfs10), aes(Var1,Var2, fill=value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = rev(terrain.colors(30)), limits = c(0,3.2))+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  #xlab("Southern population") + ylab("Northern population")+
  xlab("") + ylab("")+
  ggtitle("H. boulengeri - Central")+
  #labs(fill = "Alleles\n(log10)")+
  theme(legend.position = "none") +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 10)); HboulM

#HboulN
SFS <- readLines(paste("./",filez[4],"/",filez[4],"_MSFS.obs", sep=""))
sfs_num <- as.numeric(unlist(strsplit(SFS[2], "\t")))
sfs <- matrix(as.numeric(unlist(strsplit(SFS[3], "\t"))), nrow = sfs_num[2]+1, ncol = sfs_num[2]+1, byrow = T)
sfs10 <- log10(sfs)
sfs10[which(sfs10 == -Inf)] <- 0
rownames(sfs10) <- seq(0,sfs_num[2],1)
colnames(sfs10) <- seq(0,sfs_num[2],1)

HboulN <- ggplot(melt(sfs10), aes(Var1,Var2, fill=value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = rev(terrain.colors(30)), limits = c(0,3.2))+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  #xlab("Southern population") + ylab("Northern population")+
  xlab("") + ylab("")+
  ggtitle("H. boulengeri - North")+
  #labs(fill = "Alleles\n(log10)")+
  theme(legend.position = "none") +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 10)); HboulN
                
#HboulS
SFS <- readLines(paste("./",filez[5],"/",filez[5],"_MSFS.obs", sep=""))
sfs_num <- as.numeric(unlist(strsplit(SFS[2], "\t")))
sfs <- matrix(as.numeric(unlist(strsplit(SFS[3], "\t"))), nrow = sfs_num[2]+1, ncol = sfs_num[2]+1, byrow = T)
sfs10 <- log10(sfs)
sfs10[which(sfs10 == -Inf)] <- 0
rownames(sfs10) <- seq(0,sfs_num[2],1)
colnames(sfs10) <- seq(0,sfs_num[2],1)

HboulS <- ggplot(melt(sfs10), aes(Var1,Var2, fill=value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = rev(terrain.colors(30)), limits = c(0,3.2))+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  #xlab("Southern population") + ylab("Northern population")+
  xlab("") + ylab("")+
  ggtitle("H. boulengeri - South")+
  #labs(fill = "Alleles\n(log10)")+
  theme(legend.position = "none") +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 10)); HboulS

#HolM
SFS <- readLines(paste("./",filez[6],"/",filez[6],"_MSFS.obs", sep=""))
sfs_num <- as.numeric(unlist(strsplit(SFS[2], "\t")))
sfs <- matrix(as.numeric(unlist(strsplit(SFS[3], "\t"))), nrow = sfs_num[2]+1, ncol = sfs_num[2]+1, byrow = T)
sfs10 <- log10(sfs)
sfs10[which(sfs10 == -Inf)] <- 0
rownames(sfs10) <- seq(0,sfs_num[2],1)
colnames(sfs10) <- seq(0,sfs_num[2],1)

HolM <- ggplot(melt(sfs10), aes(Var1,Var2, fill=value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = rev(terrain.colors(30)), limits = c(0,3.2))+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  #xlab("Southern population") + ylab("Northern population")+
  xlab("") + ylab("")+
  ggtitle("Hollandichthys - Central")+
  #labs(fill = "Alleles\n(log10)")+
  theme(legend.position = "none") +
  #theme(legend.position = c(0.9, 0.65))+ 
  guides(fill = guide_colourbar(barwidth = 1, barheight = 10)); HolM

#HolS
SFS <- readLines(paste("./",filez[7],"/",filez[7],"_MSFS.obs", sep=""))
sfs_num <- as.numeric(unlist(strsplit(SFS[2], "\t")))
sfs <- matrix(as.numeric(unlist(strsplit(SFS[3], "\t"))), nrow = sfs_num[2]+1, ncol = sfs_num[2]+1, byrow = T)
sfs10 <- log10(sfs)
sfs10[which(sfs10 == -Inf)] <- 0
rownames(sfs10) <- seq(0,sfs_num[2],1)
colnames(sfs10) <- seq(0,sfs_num[2],1)

HolS <- ggplot(melt(sfs10), aes(Var1,Var2, fill=value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = rev(terrain.colors(30)), limits = c(0,3.2))+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  #xlab("Southern population") + ylab("Northern population")+
  xlab("") + ylab("")+
  ggtitle("Hollandichthys - South")+
  #labs(fill = "Alleles\n(log10)")+
  theme(legend.position = c(0.9, 0.65))+ 
  theme(legend.position = "none") +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 10)); HolS
                
#MimaM
SFS <- readLines(paste("./",filez[8],"/",filez[8],"_MSFS.obs", sep=""))
sfs_num <- as.numeric(unlist(strsplit(SFS[2], "\t")))
sfs <- matrix(as.numeric(unlist(strsplit(SFS[3], "\t"))), nrow = sfs_num[2]+1, ncol = sfs_num[2]+1, byrow = T)
sfs10 <- log10(sfs)
sfs10[which(sfs10 == -Inf)] <- 0
rownames(sfs10) <- seq(0,sfs_num[2],1)
colnames(sfs10) <- seq(0,sfs_num[2],1)

MimaM <- ggplot(melt(sfs10), aes(Var1,Var2, fill=value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = rev(terrain.colors(30)), limits = c(0,3.2))+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  #xlab("Southern population") + ylab("Northern population")+
  xlab("") + ylab("")+
  ggtitle("M. microlepis - Central")+
  #labs(fill = "Alleles\n(log10)")+
  #theme(legend.position = c(0.9, 0.65))+ 
  theme(legend.position = "none") +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 10)); MimaM

#MimaN
SFS <- readLines(paste("./",filez[9],"/",filez[9],"_MSFS.obs", sep=""))
sfs_num <- as.numeric(unlist(strsplit(SFS[2], "\t")))
sfs <- matrix(as.numeric(unlist(strsplit(SFS[3], "\t"))), nrow = sfs_num[2]+1, ncol = sfs_num[2]+1, byrow = T)
sfs10 <- log10(sfs)
sfs10[which(sfs10 == -Inf)] <- 0
rownames(sfs10) <- seq(0,sfs_num[2],1)
colnames(sfs10) <- seq(0,sfs_num[2],1)

MimaN <- ggplot(melt(sfs10), aes(Var1,Var2, fill=value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = rev(terrain.colors(30)), limits = c(0,3.2))+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  #xlab("Southern population") + ylab("Northern population")+
  xlab("") + ylab("")+
  ggtitle("M. microlepis - North")+
  #labs(fill = "Alleles\n(log10)")+
  #theme(legend.position = c(0.9, 0.65))+ 
  theme(legend.position = "none") +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 10)); MimaN

#MimaS
SFS <- readLines(paste("./",filez[10],"/",filez[10],"_MSFS.obs", sep=""))
sfs_num <- as.numeric(unlist(strsplit(SFS[2], "\t")))
sfs <- matrix(as.numeric(unlist(strsplit(SFS[3], "\t"))), nrow = sfs_num[2]+1, ncol = sfs_num[2]+1, byrow = T)
sfs10 <- log10(sfs)
sfs10[which(sfs10 == -Inf)] <- 0
rownames(sfs10) <- seq(0,sfs_num[2],1)
colnames(sfs10) <- seq(0,sfs_num[2],1)

MimaS <- ggplot(melt(sfs10), aes(Var1,Var2, fill=value)) + 
  geom_raster() + 
  scale_fill_gradientn(colours = rev(terrain.colors(30)), limits = c(0,3.2))+
  scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
  #xlab("Southern population") + ylab("Northern population")+
  xlab("") + ylab("")+
  ggtitle("M. microlepis - South")+
  #labs(fill = "Alleles\n(log10)")+
  #theme(legend.position = c(0.9, 0.65))+ 
  theme(legend.position = "none") +
  guides(fill = guide_colourbar(barwidth = 1, barheight = 10)); MimaS

library(grid)
library(gridExtra)
library(lattice)
t <- textGrob("")
grid.arrange(MimaS, MimaM, MimaN, 
             HboulS, HboulM, HboulN,
             HolS, HolM, t,
             BrycoS, BrycoM, t,
             nrow = 4)

