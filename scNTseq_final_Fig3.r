dat.fig3b <- readRDS("Fig3b.rds")

library("RColorBrewer")
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

options(repr.plot.width=10,repr.plot.height=2.5)
p1 <- ggplot(dat.fig3b,aes(umap_0,umap_1)) + geom_point(aes(col=time),size=.3,alpha=.5) + 
     scale_color_manual(values = brewer.pal(n=9,name="Set1")[1:5]) + xlab("UMAP1") + ylab("UMAP2") +
     guides(colour = guide_legend(override.aes = list(size=2,alpha=1)))
p2 <- ggplot(dat.fig3b,aes(umap_0,umap_1)) + geom_point(aes(col=early),size=.3,alpha=.5) + 
     scale_color_gradientn(colours = brewer.pal(n=9,name="YlOrRd")) + xlab("UMAP1") + ylab("UMAP2")
p3 <- ggplot(dat.fig3b,aes(umap_0,umap_1)) + geom_point(aes(col=late),size=.1,alpha=.5) + 
     scale_color_gradientn(colours = brewer.pal(n=9,name="YlOrRd")) + xlab("UMAP1") + ylab("UMAP2")
plot_grid(p1,p2,p3,nrow = 1)

dat.fig3c <- readRDS("Fig3c.rds")

options(repr.plot.width=4,repr.plot.height=4.5)
ggplot(dat.fig3c,aes(variable,gene)) + geom_point(aes(size=logP,col=significant),alpha=.8) + 
 scale_colour_manual(values=c("red","grey","blue")) + 
  xlab("") + ylab("") + theme_bw()

dat.fig3d <- readRDS("Fig3d.rds")

options(repr.plot.width=15,repr.plot.height=2.5)
p1 <- ggplot(dat.fig3d,aes(umap_0,umap_1)) + geom_point(aes(col=`Jun (22g)`),size=.3,alpha=.5) + 
     scale_color_gradientn(colours = rev(brewer.pal(n=11,name="Spectral"))) + xlab("UMAP1") + ylab("UMAP2") 

p2 <- ggplot(dat.fig3d,aes(umap_0,umap_1)) + geom_point(aes(col=`Mef2d (109g)`),size=.3,alpha=.5) + 
     scale_color_gradientn(colours = rev(brewer.pal(n=11,name="Spectral"))) + xlab("UMAP1") + ylab("UMAP2")

p3 <- ggplot(dat.fig3d,aes(umap_0,umap_1)) + geom_point(aes(col=`Maff_extended (848g)`),size=.3,alpha=.5) + 
     scale_color_gradientn(colours = rev(brewer.pal(n=11,name="Spectral"))) + xlab("UMAP1") + ylab("UMAP2")

plot_grid(p1,p2,p3,nrow = 1,align="v")


