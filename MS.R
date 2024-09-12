rm(list=ls())
gc()
setwd('/home/shpcv2_kvcka/data/20230426_PDAC/Macro_Stromal')
getwd()
MS <- readRDS('Macro_Stromal_raw.rds') 

library(harmony)
MS@active.assay
MS_harmony <- SCTransform(object = MS,vars.to.regress = "percent.MT")

system.time({MS_harmony <- RunHarmony(MS_harmony, group.by.vars = "orig.ident")})
MS <- MS_harmony
rm(MS_harmony)
DimPlot(MS,reduction = 'harmony',group.by = 'orig.ident')
ElbowPlot(MS, ndims = 50)


MS@graphs$graphname
MS <- RunTSNE(MS, reduction = "harmony", dims = 1:40)

MS <- FindNeighbors(MS, reduction = "harmony", dims = 1:40)
MS <- FindClusters(MS,resolution = 0.8,graph.name = "RNA_snn")

pdf('20230506_MS_tsne.pdf')
DimPlot(MS,reduction = 'tsne',label=T)
dev.off()
table(MS$RNA_snn_res.0.8)
FeaturePlot(MS,features = 'SLC12A5',order=T,reduction = 'tsne',pt.size = 0.5,label=T)
FeaturePlot(MS,features = 'CD163',order=T,reduction = 'tsne',pt.size = 1)
mysub <- FindSubCluster(
  MS,
  graph.name='RNA_snn', 
  cluster='4', 
  subcluster.name = "sub.cluster3_MS",
  resolution = 0.5,
  algorithm = 1                           
)
Idents(mysub) <- 'sub.cluster3_MS'
pdf('20230506_sub.cluster3_MS_tsne.pdf')
DimPlot(mysub,reduction = 'tsne',label=T)
dev.off()
table(mysub@active.ident)
FeaturePlot(mysub,features = 'SLC12A5',order=T,reduction = 'tsne',pt.size = 0.5,label=T)
MS <- mysub
rm(mysub)
save(MS,file='20230507_MS.Rdata')
unique(MS@active.ident)
paste0(sprintf("'%s'", unique(MS@active.ident)), collapse = ",")
FeaturePlot(MS,features = c('PECAM1','CDH5','CD34','VWF','KDR','PLVAP'),reduction = 'tsne',pt.size = 0.7,label=T)
new.cluster.ids  <- c('mM2-TAMs','M1-TAMs','bM2-TAMs','mM2-TAMs','M1-TAMs',
                      'bM2-TAMs','Endothelial_cells','mM2-TAMs','bM2-TAMs','myCAFs',
                      'mM2-TAMs','Stellate_cells','mM2-TAMs','Acinar_cells','Acinar_cells',
                      'mM2-TAMs','14','iCAFs','myCAFs','mM2-TAMs',
                      'Acinar_cells','Acinar_cells','bM2-TAMs','mM2-TAMs','19')
table(MS@active.ident)
names(new.cluster.ids) <- levels(MS)
MS <- RenameIdents(MS, new.cluster.ids)
table(MS@active.ident)
MS$SLC12A5_MS <- MS@active.ident
MS <- subset(MS,subset = SLC12A5_MS!='14')
MS <- subset(MS,subset = SLC12A5_MS!='19')
pdf('20230507_MS_annotated.pdf')
DimPlot(MS,reduction = 'tsne',label = T)
dev.off()
save(MS,file='20230507_MS_annotated.Rdata')


mydata <- MS@reductions$tsne@cell.embeddings 
mymeta <- as.data.frame(MS@meta.data)   
mydata <- cbind(mydata,mymeta)
head(mydata)
allcolour=c("#20B2AA","#0000FF","#DC143C",'#D3D3D3','#9400B3')
allcolour=c("#E64B35FF","#4DBBD5FF","#9370DB","#1E90FF","#FFA500","#228B22","#5F9EA0",'#4682B4',"#800080","#A0522D",'#87ceeb',"#D2691E","#87CEEB","#40E0D0","#5F9EA0",   "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
p <- ggplot(mydata,aes(x= tSNE_1 , y = tSNE_2 ,color = SLC12A5_MS)) +  geom_point(size = 0.5 , alpha =1 )  +  scale_color_manual(values = allcolour)+theme(legend.position = 'none')
p1 <- p+theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.border = element_blank(), 
              axis.title = element_blank(),  
              axis.text = element_blank(), 
              axis.ticks = element_blank(),
              panel.background = element_rect(fill = 'white'), 
              plot.background=element_rect(fill="white"))
p2 <- p1+theme(
  legend.title = element_blank(), 
  legend.key=element_rect(fill='white'), 
  legend.text = element_text(size=20), 
  legend.key.size=unit(1,'cm') ) +  
  guides(color = guide_legend(override.aes = list(size=5))) 
p3 <- p2 + 
  geom_segment(aes(x = min(mydata$tSNE_1) , y = min(mydata$tSNE_2) ,
                   xend = min(mydata$tSNE_1) +8, yend = min(mydata$tSNE_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(mydata$tSNE_1)  , y = min(mydata$tSNE_2)  ,  
                   xend = min(mydata$tSNE_1) , yend = min(mydata$tSNE_2) + 8),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(mydata$tSNE_1) +3, y = min(mydata$tSNE_2) -2.5, label = "tSNE_1",
           color="black",size = 3, fontface="bold" ) +  
  annotate("text", x = min(mydata$tSNE_1) -2.5, y = min(mydata$tSNE_2) + 3, label = "tSNE_2",
           color="black",size = 3, fontface="bold" ,angle=90) 

library(tidyverse)
cell_type_med <- mydata %>%
  group_by(SLC12A5_MS) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )
library(ggrepel)

p4 <- p3 +
  geom_label_repel(aes(label=SLC12A5_MS), fontface="bold",data = cell_type_med,
                   point.padding=unit(0.5, "lines"),size=5) +
  theme(legend.position = "none")
#ggsave(("20230110_immu_dimplot.pdf"),p4,height=18,width=22,unit="cm")    

pdf('20230507_MS_SLC12A5_MS_annotated.pdf',height =6,width = 6)
p4
dev.off()




pdf('20230510_VCAM_pathway.pdf',height =3,width =3)
VlnPlot(MS,features = c('VCAM1'),assay = 'RNA',ncol = 1,pt.size = 0)+scale_color_npg()+ theme(text = element_text(size =10))
dev.off()


library(Seurat)
library(ggplot2)
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45), axis.ticks.x = element_line())
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}
pdf('20230510_VCAM_pathway.pdf',height =5,width =3)
StackedVlnPlot(MS, c('VCAM1','ITGA4','ITGA9','ITGB1'), pt.size=0, cols=my36colors,assay='RNA')
dev.off()

my36colors <- c('#53A85F', '#F1BB72', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
