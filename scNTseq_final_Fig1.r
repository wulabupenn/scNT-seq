library(ggplot2)
library(cowplot)
library(dplyr)
library(tidyr)

dat.fig1b <- readRDS("Fig1b.rds") 

head(dat.fig1b,n=2)

mix3 <- dat.fig1b %>% group_by(species,method) %>% dplyr::summarise(n=n())
mix3.1 <- mix3 %>% filter(species=="human")
mix3.2 <- mix3 %>% filter(species=="mouse")
mix3.3 <- mix3 %>% filter(species=="mixed")

options(repr.plot.width=9,repr.plot.height=3)
ggplot(dat.fig1b,aes(human,mouse)) + geom_point(aes(color=species),alpha=0.8,size=0.7,shape=1) + 
  scale_colour_brewer(palette = "Set1") +
  xlab("K562") + ylab("mESC")+
  geom_text(data=mix3.2, aes(label=paste(species,n,sep=" ")), 
            x=Inf, y=Inf, hjust=1.1, vjust=1.2,size=3,
            colour=RColorBrewer::brewer.pal(8,"Set1")[1], inherit.aes=FALSE, parse=FALSE)+
  geom_text(data=mix3.1, aes(label=paste(species,n,sep=" ")), 
            x=Inf, y=Inf, hjust=1.1, vjust=2.7,size=3,
            colour=RColorBrewer::brewer.pal(8,"Set1")[2], inherit.aes=FALSE, parse=FALSE)+
  geom_text(data=mix3.3, aes(label=paste(species,n,sep=" ")), 
            x=Inf, y=Inf, hjust=1.1, vjust=4.2,size=3,
            colour=RColorBrewer::brewer.pal(8,"Set1")[3], inherit.aes=FALSE, parse=FALSE)+
  
  facet_wrap(~method,scales="free",nrow=1) + 
  theme_bw() + 
  theme(legend.position="none",legend.title=element_blank(),panel.grid.minor = element_blank())

dat.fig1c <- readRDS("Fig1c.rds")

options(repr.plot.height=4,repr.plot.width=8)
ggplot(dat.fig1c,aes(type,percent2,fill=sample)) + 
 geom_bar(stat = "identity", position="dodge",width=.7) + 
 ylab("percentage (%)") + xlab("") + 
 geom_text(aes(label=round(percent2,2)), position=position_dodge(width=0.7), vjust=-0.25,size=2.3)  + 
 scale_fill_brewer("",palette = "Set1",labels=c("Untreated","TFEA/NaIO4"))+ 
 theme_bw() + ggtitle("K562")

dat.fig1d <- readRDS("Fig1d.rds")

options(repr.plot.width=2,repr.plot.height=3)
ggplot(dat.fig1d, aes(type,label.UMI.percent))  + geom_boxplot() + xlab("")

dat.fig1e <- readRDS("Fig1e.rds")

gene.plot.p1 <- dat.fig1e %>% filter(defined.base=="T")
gene.plot.p2 <- dat.fig1e %>% filter(defined.base=="C")

p.actg1 <- ggplot(dat.fig1e,aes(bp,index)) + 
   geom_point(data=gene.plot.p1,aes(bp,index),shape=1,col="grey60",alpha=.5,size=0.3) + 
   scale_color_gradient("T>C coverage",low = "blue",high = "red") + 
   geom_point(data=gene.plot.p2,aes(bp,index,col=C),shape=4,size=0.5) + scale_y_discrete(labels=NULL) + 
    theme_bw() + ggtitle("ACTG1") + ylab("41 UMIs (19 label, 22 unlabel)")

read.plot <- readRDS("Fig1e_bottom.rds")

read.plot.p1 <- read.plot %>% filter(base=="T") %>% droplevels
read.plot.p2 <- read.plot %>% filter(base=="C") %>%  droplevels

p.example <- ggplot(read.plot,aes(bp,read)) + 
   geom_point(data=read.plot.p1,aes(bp,read),shape=1,col="grey60",size=2) + geom_line(linetype="dashed",size=.1) +
   geom_point(data=read.plot.p2,aes(bp,read,col=base),shape=4,size=3,col="blue") + scale_y_discrete(labels=NULL)  + theme_bw() +
   ggtitle("Reads for 2nd UMI")

options(repr.plot.height=3.5,repr.plot.width=10)
plot_grid(p.actg1,p.example,nrow=1,rel_widths = c(1.5,1))


