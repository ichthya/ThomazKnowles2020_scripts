library(ggplot2)
#library(cowplot)

setwd('~/Dropbox/Comparative/Fastsimcoal/')

tdiv <- read.csv('Fatsimcoal_TDIV_plot_Dec2019.csv')[1:10,]
tdiv$Break_f <- factor(tdiv$Break, levels = c("North", "Central", "South"))
tdiv$Species_f <- factor(tdiv$Species, 
                         levels = c("Bryconamericus",
                                    "Hollandichthys",
                                    "H. boulengeri",
                                    "M. microlepis"))
sea <- read.csv('sea_level_Pleistocene.csv')
sea$sea <- "Sea Level (m)"

max_sea <- max(tdiv$max_TDIV)
sea$kya <- sea$Age_Ma*10^6
sea_comp <- subset(sea, kya < max_sea)

col <- c("Central" = "black", "North" = "steelblue", "South" = "red4")

p <- ggplot(tdiv, aes(x=tdiv$Species_f,
                      y = TDIV, ymin=min_TDIV, ymax=max_TDIV, 
                      colour=factor(Break_f)))+
     geom_pointrange(size = 0.8)+
     #geom_pointrange(data = tdiv, aes(x = tdiv$comb, y = TADM, ymin = min_TADM, ymax = max_TADM), size = 0.4)+
     #geom_ribbon(aes(ymin = min_MIG, max = max_MIG), linetype = "dotted", size = 0.4)+
     coord_flip()+
     scale_colour_manual(values = col, guide = F)+
     theme_bw()+
     theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
           axis.text.y = element_text(face = "italic", size = 12),
           plot.margin = margin(t=5,r=5,b=0,l=7),
           strip.text.y = element_text(size = 12, face = "bold"),
           axis.ticks.x=element_blank())+
     expand_limits(y= c(0, max_sea))+
     scale_y_continuous(expand=c(0,0), breaks = c(0, 50000,100000,150000,200000,250000))+
     facet_grid(Break_f ~ ., scales = "free_y", space = "free_y")+
     labs(y = "Time (Kya)", x=""); p

s <- ggplot(sea_comp, aes(x=kya, y=Sealevel, ymin = -125, ymax = 15, xmax = max_sea))+
     geom_line(color="darkblue", lwd = 2)+
     geom_ribbon(aes(ymin = -125, ymax= Sealevel), fill="lightblue")+
     scale_x_continuous(expand=c(0,0), breaks = c(0, 50000,100000,150000,200000,250000), 
                        labels = c(0,50,100,150,200,250))+
     scale_y_continuous(expand=c(0,0))+
     theme_bw()+
     theme(axis.title.y=element_blank(), plot.margin = margin(t=3,r=5,b=5,l=80),
           axis.text = element_text(size = 12),
           #axis.text.x = element_text(angle = 45),
           strip.text.y = element_text(size = 12, face = "bold"))+
     facet_grid(sea ~ .)+
     labs(x = "Time (Kya)", y="Sea level"); s

library(cowplot)
plot_grid(p, s, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2,
          rel_heights = c(2.3,1))

