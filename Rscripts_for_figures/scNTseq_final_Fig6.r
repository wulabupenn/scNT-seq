dat.fig6b <- readRDS("Fig6b.rds")

library(mgcv)

library(dplyr)
library(tidyr)

sampleID <- c('K562-4SU','K562-NTC-TFEA','K562-4SU-TFEA','K562-4SU-N9','K562-NTC-TFEA-N9','K562-4SU-TFEA-N9')

predict.nGene <- sapply(sampleID,function(u){
    model1 <- loess(data = dat.fig6b %>% filter(orig.ident == u),formula = nGene ~ reads)
    tmp <- predict(model1,50000)
    return(round(as.numeric(tmp),0))
})
df.extra <- data.frame(reads=rep(50000,6),nGene=predict.nGene)

df.extra$orig.ident <- rownames(df.extra)
library(RColorBrewer)
my.col <- brewer.pal("Set1",n=6)[c(1,3,2,4,6,5)]
dat.fig6b$orig.ident <- factor(dat.fig6b$orig.ident,
                             levels=c('K562-4SU','K562-NTC-TFEA','K562-4SU-TFEA',
                                      'K562-4SU-N9','K562-NTC-TFEA-N9','K562-4SU-TFEA-N9'))

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

options(repr.plot.width=6,repr.plot.height=4)
ggplot(dat.fig6b,aes(reads,nGene,col=orig.ident)) + geom_point(size=.2,shape=19,alpha=.3) + 
  geom_smooth(method="loess",se=F,size=.5,linetype=2,aes(group=orig.ident),col="black") +
  scale_color_manual("Method",values=my.col,labels=c("4SU","TFEA","4SU-TFEA","4SU-N9","TFEA-N9","4SU-TFEA-N9")) +
  scale_x_continuous(breaks=seq(0,80000,25000),labels=seq(0,80000,25000),limits=c(0,80000)) + ylim(0,9000) + 
  geom_vline(xintercept = 50000,col="grey30",linetype="dashed") + 
  geom_point(data = df.extra, col = my.col) +
  geom_text(data=df.extra,aes(label=nGene),col="black",hjust=-.1) +
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  guides(colour = guide_legend(ncol=1,override.aes=list(size=3,alpha=1)))

predict.nUMI <- sapply(sampleID,function(u){
    model1 <- lm(data = dat.fig6b %>% filter(orig.ident == u),formula = nUMI ~ reads)
    tmp <- predict(model1,list(reads=seq(30000,50000,10000)))
    return(round(as.numeric(tmp),0))
})

df.extra2 <- data.frame(reads=rep(50000,6),nUMI=predict.nUMI[3,])
df.extra2$orig.ident <- rownames(df.extra2)

options(repr.plot.width=6,repr.plot.height=4)
ggplot(dat.fig6b,aes(reads,nUMI,col=orig.ident)) + geom_point(size=.3,shape=19,alpha=.5) + 
  geom_smooth(method="lm",formula = y ~ x,se=F,size=.5,linetype=2,aes(group=orig.ident),col="black") +
  scale_color_manual("Method",values=my.col,labels=c("4SU","TFEA","4SU-TFEA","4SU-SW","TFEA-SW","4SU-TFEA-SW")) +
  scale_x_continuous(breaks=seq(0,80000,25000),labels=seq(0,80000,25000),limits=c(0,80000)) + ylim(0,50000) +
  geom_vline(xintercept = 50000,col="grey30",linetype="dashed") + 
  geom_point(data = df.extra2, col = my.col) +
  geom_text(data=df.extra2,aes(label=nUMI),col="black",hjust=-.1) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  guides(colour = guide_legend(ncol=1,override.aes=list(size=3,alpha=1)))

dat.fig6c <- readRDS("Fig6c.rds")

library(MASS)
library(ggplot2)
library(viridis)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# new
dat <- dat.fig6c$new %>% dplyr::rename(x=WO,y=N9)
dat$density <- get_density(dat$x, dat$y, n = 1000)
p.new <- ggplot(dat) + geom_point(aes(x, y, color = density),size=.3) + 
 scale_color_distiller(palette = "Spectral") +
 xlab("4SU/TFEA") + ylab("4SU/TFEA/2nd SS") +
 annotate("text",x=0.2,y=.9,label="0.91",size=4,color="red") 

# old
dat <- dat.fig6c$old  %>% dplyr::rename(x=WO,y=N9)
dat$density <- get_density(dat$x, dat$y, n = 1000)
p.old <- ggplot(dat) + geom_point(aes(x, y, color = density),size=.3) + scale_color_distiller(palette = "Spectral") +
 xlab("4SU/TFEA") + ylab("4SU/TFEA/2nd SS") +
 annotate("text",x=0.4,y=4.9,label="0.92",size=4,color="red") + xlim(0,6) + ylim(0,6)

# ntr
dat <- dat.fig6c$ntr %>% dplyr::rename(x=WO,y=N9)
dat$density <- get_density(dat$x, dat$y, n = 1000)
p.ntr <- ggplot(dat) + geom_point(aes(x, y, color = density),size=.3) + scale_color_distiller(palette = "Spectral") +
 xlab("4SU/TFEA") + ylab("4SU/TFEA/2nd SS") +
 annotate("text",x=0.1,y=0.9,label="0.88",size=4,color="red") 

options(repr.plot.width=11,repr.plot.height=3)
plot_grid(p.new,p.old,p.ntr,nrow=1)

dat.fig6d <- readRDS("Fig6d.rds")

options(repr.plot.width=7,repr.plot.height=6)
p1 <- ggplot(dat.fig6d,aes(PC_1,PC_2,col=Phase)) + geom_point(size=.4,alpha=.4)+ theme_cowplot() + 
  guides(colour = guide_legend(title="Cluster",override.aes = list(size=4))) + xlab("PC 1") + ylab("PC 2")
p2 <- ggplot(dat.fig6d,aes(PC_1,PC_2,col=orig.ident)) + geom_point(size=.4,alpha=.8) + theme_cowplot() +
  scale_color_manual(values=c('#377EB8','#984EA3','#FF7F00'),
                    labels=c("4sU/TFEA","4sU/2nd SS","4sU/TFEA/2nd SS"))+
  guides(colour = guide_legend(title="Sample",override.aes = list(size=4)))  + xlab("PC 1") + ylab("PC 2")
plot_grid(p1,p2,align="v",nrow=2)

dat.6e <- readRDS("Fig6e.rds")

options(repr.plot.height=4,repr.plot.width=12)
ggplot(dat.6e,aes(orig.ident,value,fill=orig.ident))+
 geom_violin(trim=T,scale="width") + 
 geom_boxplot(width=0.1) +
 facet_grid(.~ Var1,scales="free_y") + 
 theme(strip.background = element_blank(), strip.placement = "outside",#axis.text.y=element_blank(), axis.title.y=element_blank(),
      strip.text.y = element_text(colour = "red", angle = 360,size=10),legend.position="none",
      panel.grid=element_blank(), panel.border=element_blank()) +
     theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1,size=rel(0.9))) + xlab("") + ylab("ratio") + 
     scale_fill_manual(values=c('#377EB8','#FF7F00'),labels=c("4sU/TFEA","4sU/TFEA/2nd SS")) 

dat.6f <- readRDS("Fig6f.rds")

options(repr.plot.height=3,repr.plot.width=7)
ggplot(dat.6f,aes(orig.ident,value,fill=transcript))+
       geom_violin(scale="width",trim=T,alpha=0.8,adjust=1,width=.7)+
       facet_grid(gene~ Phase,scales="free_y") + 
       theme(strip.background = element_blank(), strip.placement = "outside",
             strip.text.y = element_text(colour = "red", angle = 360,size=10),legend.position="right",
             panel.grid=element_blank(), panel.border=element_blank()) +
       theme(axis.text.x = element_text(angle = 45,hjust=1,vjust=1,size=rel(0.9))) + xlab("") + ylab("normalized expression level") + 
       scale_fill_manual(values=rev(brewer.pal(3, "Set1"))) 


