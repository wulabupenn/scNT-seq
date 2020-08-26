dat.fig2b <- readRDS("Fig2b.rds")

library("Seurat") # here seurat version #2.3.4

packageVersion("Seurat") # 2.3.4

library("RColorBrewer")

options(repr.plot.height=4,repr.plot.width=4)
DimPlot(object = dat.fig2b,reduction.use = "umap",do.label =T, no.legend = T, do.return = TRUE, pt.size = 0.1,
              cols.use = brewer.pal("Set1",n=9))

options(repr.plot.height=3.2,repr.plot.width=4)
DimPlot(object = dat.fig2b,reduction.use = "umap",do.label =F, no.legend = F, do.return = TRUE, pt.size = 0.1,group.by="time",
                 cols.use = brewer.pal("Spectral",n=5))

dat.fig2c <- readRDS(file="Fig2c.rds")

library("Seurat")
library(ggplot2)
library(RColorBrewer)
library(cowplot)
theme_set(theme_cowplot())

# switch to Seurat 3 here
packageVersion("Seurat") # 3.1.4

plot_time_cluster <- function(mat) {
     m.t <- CreateSeuratObject(counts = mat,min.cells = 10)
     m.t@meta.data <- dat.fig2c$cell.info[rownames(m.t@meta.data),]
     m.t@meta.data <- droplevels(m.t@meta.data)
     m.t <- FindVariableFeatures(m.t, selection.method = "vst", nfeatures = 2000,verbose = F)
     m.t <- ScaleData(m.t,verbose = F)
     m.t <- RunPCA(m.t,npcs = 20,features =  VariableFeatures(object = m.t),verbose=FALSE)
     dff <- cbind(m.t@reductions$pca@cell.embeddings,m.t@meta.data)
     p1 <- ggplot(dff,aes(PC_1,PC_2,col=as.factor(time))) + geom_point(size=.5,alpha=.3) +
           scale_color_manual(values=brewer.pal("Set1",n=3),"Time") + 
           geom_density_2d(size=.05) +
           xlab("PC 1") + ylab("PC 2") + 
           facet_wrap(~cluster,scales="free")  +
           theme(text = element_text(size=12),strip.background = element_blank(),
           strip.placement = "outside")
  return(p1)
}

options(repr.plot.height=5,repr.plot.width=11)
p1 <- plot_time_cluster(dat.fig2c$new) + ggtitle("new")
p2 <- plot_time_cluster(dat.fig2c$old)+ ggtitle("old")
p3 <- plot_time_cluster(dat.fig2c$total)+ ggtitle("total")
p4 <- plot_time_cluster(dat.fig2c$ntr)+ ggtitle("NTR")

plot_grid(p1,p2,p3,p4,nrow=2)

dat.fig2d <- readRDS("Fig2d.rds")

library(RColorBrewer)
options(repr.plot.height=6,repr.plot.width=5)
ggplot(dat.fig2d,aes(time,value,col=cluster.name)) + 
 geom_point() + 
 geom_line() + 
 facet_grid(Var1~type,scales = "free_y") + 
 xlab("time (min)") + ylab("normalized expression") + scale_color_manual(values=brewer.pal("Set1",n=9)) + 
 scale_x_continuous(breaks=c(0,15,30,60,120)) +
 theme(axis.text.x = element_text(angle=45,hjust=1))

dat.fig2e <- readRDS("Fig2e.rds")

options(repr.plot.height=4,repr.plot.width=14)
pheatmap::pheatmap(dat.fig2e$induced,show_colnames = T,cluster_rows = T,cluster_cols = F,
                   color = rev(colorRampPalette(brewer.pal(9,"Spectral"))(50)),
                  fontsize = 7,clustering_method ="ward.D",border_color ="NA",
                  gaps_col = c(45))

options(repr.plot.height=10,repr.plot.width=14)
pheatmap::pheatmap(dat.fig2e$noninduced,show_colnames = T,cluster_rows = T,cluster_cols = F,
                   color = rev(colorRampPalette(brewer.pal(9,"Spectral"))(50)),
                  fontsize = 7,clustering_method ="ward.D",border_color ="NA",
                  gaps_col = c(45))

dat.fig2f <- readRDS("Fig2f.rds")

options(repr.plot.width=10,repr.plot.height=5)
ggplot(dat.fig2f,aes(cluster_time,AUC,fill = cluster2)) + 
       geom_boxplot(outlier.size=.1) + 
       facet_grid(TF~type,scales = "free_y") +
       xlab("") + 
       ylab("AUC value") + 
       theme(axis.text.x=element_text(angle=45,hjust=1),strip.background = element_blank(),
             strip.text.y = element_text(size = 8, colour = "red"),
            legend.position = "none")


