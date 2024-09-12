setwd("/home/shpcv2_kvcka/data/20230426_PDAC/Immune")
getwd()

Immune <- obj.combined
rm(obj.combined)
gc()

Immune@graphs$graphname
Immune <- RunTSNE(Immune, reduction = "harmony", dims = 1:40)
Immune <- FindNeighbors(Immune, reduction = "harmony", dims = 1:40)
Immune <- FindClusters(Immune,resolution = 0.8,graph.name = "RNA_snn")
DimPlot(Immune,reduction = 'tsne')
save(Immune,file='20230428_Immune.Rdata')


setwd("/home/shpcv2_kvcka/data/20230426_PDAC/Immune")
load('20230428_Immune.Rdata')
pdf('tsne.pdf', width = 10,height = 10)
DimPlot(Immune,reduction = 'tsne',label = T)
dev.off()

T_cell <- c('CD3D','IL7R','CD3E')
CD4_T <- c('IL7R','LTB','LDHB','CD4')
CD8_T <- c('NKG7','GNLY','CCL5','CD8A')
NK <- c('FCGR3A','NKG7','GNLY','GZMB','KLRB1')
Myeloid_cell <- c('CD68','LYZ','CD14','IL3RA','LAMP3','CLEC4C','TPSAB1')
Macro <- c('CD14','CD68','CD163')
B_cell <- c('CD79A','MS4A1','CD79B','CD19')
DC <- c('TSPAN13','GPR183','CD1C','GSN')
Mast <- c('CPA3','TPSAB1','TPSB2')
Plasma <- c("IGJ",'CD79A','IGKC','IGHG3','IGHG1')
Treg <- c('FOXP3','CTLA4','IL10','TNFRSF18','CCR8','IKZF4','IKZF2','TNFRSF4','IL2RA')
Neutrophils <- c('FPR1','G0S2','S100A8','S100A9','FCGR3B','CXCR2','CXCR1')
gene_list = list(
  T_cell= T_cell,CD4_T=CD4_T,CD8_T=CD8_T,
  NK = NK,
  Macro=Macro,B_cell=B_cell,DC=DC,Mast=Mast,
  Plasma=Plasma,Treg=Treg,Neutrophils=Neutrophils
)
genes_to_check=unique(unlist(gene_list))

pdf('tsne_markers.pdf',width = 10,height =8)
DotPlot(Immune,assay = "RNA",cols = rep(c("darkblue", "darkred",'darkorange','black','darkgreen'),8),features = genes_to_check, split.by = "RNA_snn_res.0.8") + coord_flip()+ RotatedAxis()
dev.off()

FeaturePlot(Immune,features = c("PCNA","MKI67"),reduction = "tsne")

new.cluster.ids <- c('CD4_T',"CD8_T","Neutrophils", "M2-like TAM","M2-like TAM","M1-like TAM",
                     'B_cells',"Mast","NK","CD8_T","Plasma",
                     "DC","Treg","NK","Cycle","Neutrophils",
                     "unknown","unknown","unknown","CD8_T","unknown",
                     "unknown","unknown","M2-like TAM","M2-like TAM","unknown" ,
                     "unknown" , "unknown" , "unknown" , "unknown")
table(Immune@active.ident)
names(new.cluster.ids) <- levels(Immune)
Immune <- RenameIdents(Immune, new.cluster.ids)
table(Immune@active.ident)
DimPlot(Immune,reduction = 'tsne',label = T)
Immune$subcelltype <- Immune@active.ident



mydata <- Immune@reductions$tsne@cell.embeddings 
mymeta <- as.data.frame(Immune@meta.data)   
mydata <- cbind(mydata,mymeta)
head(mydata)
allcolour=c("#20B2AA","#0000FF","#DC143C",'#D3D3D3','#9400B3')
allcolour=c("#E64B35FF","#4DBBD5FF","#FFA500","#9370DB",'#D2B48C',"#1E90FF","#808000",'#A9A9A9',"#228B22","#7B68EE",'#4682B4',"#800080","#A0522D",'#87ceeb',"#D2691E","#87CEEB","#40E0D0","#5F9EA0",   "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
p <- ggplot(mydata,aes(x= tSNE_1 , y = tSNE_2 ,color = mbcelltype)) +  geom_point(size = 0.5 , alpha =1 )  +  scale_color_manual(values = allcolour)+theme(legend.position = 'none')
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
  group_by(mbcelltype) %>%
  summarise(
    tSNE_1 = median(tSNE_1),
    tSNE_2 = median(tSNE_2)
  )
library(ggrepel)
cell_type_med[8,"tSNE_2"] <- -40
cell_type_med[8,"tSNE_1"] <- 1.5
p4 <- p3 +
  geom_label_repel(aes(label=mbcelltype), fontface="bold",data = cell_type_med,
                   point.padding=unit(0.5, "lines"),size=5) +
  theme(legend.position = "none")
#ggsave(("20230110_immu_dimplot.pdf"),p4,height=18,width=22,unit="cm")    

pdf('20230508_Immune_mbcelltype_annotated.pdf',height=6,width=6)
p4
dev.off()

saveRDS(Immune,'20230429_Immune_annotated.rds')
rm(list=ls())

setwd("/home/shpcv2_kvcka/data/20230426_PDAC/Immune")
getwd()
Immune <- readRDS('20230429_Immune_annotated.rds')


raw.data <- as.matrix(GetAssayData(Immune,slot = 'counts'))
Immune@meta.data[['SLC12A5_subcelltype']] <- c()
for (i in 1:length(colnames(raw.data))) {
  if (raw.data['SLC12A5',colnames(raw.data)[i]]>0 & Immune@meta.data[["subcelltype"]][i] == 'M2-like TAM'){
    Immune@meta.data[["SLC12A5_subcelltype"]][i] <- "SLC12A5_high_M2"
  }
  else if(Immune@meta.data[["subcelltype"]][i] == 'M2-like TAM') {
    Immune@meta.data[["SLC12A5_subcelltype"]][i] <- "SLC12A5_low_M2"
  }
  else{
    Immune@meta.data[["SLC12A5_subcelltype"]][i] <- as.character(Immune@meta.data[["subcelltype"]][i])
  }
}

Idents(Immune) <- "SLC12A5_subcelltype"
DimPlot(Immune,reduction = 'tsne',order = T,pt.size = 1)
saveRDS(Immune,'20230430_Immune_annotated.rds')

library(ggplot2)
color <- c('lightgrey','darkred')
pdf('20230428_SLC12A5_FeaturePlot.pdf')
FeaturePlot(Immune,order=T,reduction='tsne',features = c('SLC12A5'),ncol = 1,cols=color,pt.size=1)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))
dev.off()

color <- c('lightgrey','darkred')
pdf('20230508_SIGLEC1_FeaturePlot.pdf')
FeaturePlot(Immune,order=T,reduction='tsne',features = c('SIGLEC1'),ncol = 1,cols=color,pt.size=1)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+ theme(text = element_text(size =30))
dev.off()



Idents(Immune) <- 'subcelltype'
DimPlot(Immune,label=T,reduction = 'tsne')
mysub <- FindSubCluster(
  Immune,
  graph.name='RNA_snn', 
  cluster='M2-like TAM', 
  subcluster.name = "sub.cluster2",
  resolution = 0.8,
  algorithm = 1                           
)
Idents(mysub) <- 'sub.cluster2'
pdf('20230503_M2sub.cluster2.pdf',width = 30,height = 30)
DimPlot(mysub,label=T,reduction = 'tsne',pt.size = 2)
dev.off()
pdf('20230503_M2sub.cluster2_SLC12A5.pdf',width = 30,height = 30)
FeaturePlot(mysub,features = 'SLC12A5',pt.size = 2,order = T,reduction = 'tsne',label = T)
dev.off()
Idents(mysub)
mysub@active.ident
levels(mysub)
new.cluster.ids <- c("SLC12A5_low_M2","M1-like TAM","Cycle","SLC12A5_high_M2" ,
                     "CD4_T","SLC12A5_low_M2","DC","Treg",          
                     "unknown","SLC12A5_low_M2","SLC12A5_high_M2","SLC12A5_low_M2", 
                     "CD8_T","SLC12A5_low_M2","B_cells","NK",            
                     "SLC12A5_low_M2","Neutrophils","SLC12A5_low_M2","SLC12A5_high_M2", 
                     "Plasma","Mast","SLC12A5_low_M2", "SLC12A5_low_M2",
                     "SLC12A5_low_M2")
table(mysub@active.ident)
names(new.cluster.ids) <- levels(mysub)
mysub <- RenameIdents(mysub, new.cluster.ids)
table(mysub@active.ident)
DimPlot(mysub,reduction = 'tsne',label = T)
mysub$SLC12A5_macro <-  mysub@active.ident
Immune <-  mysub
DimPlot(Immune,reduction = 'tsne',label = T)
saveRDS(Immune,'20230503_Immune_annotated_with_SLC12A5_macro.rds')

rm(list=ls())
setwd("/home/shpcv2_kvcka/data/20230426_PDAC/Immune")
getwd()
Immune <- readRDS('20230503_Immune_annotated_with_SLC12A5_macro.rds')




library(Seurat)
library(patchwork)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(paletteer)
mytheme = theme(plot.title = element_text(size = 24,color="black",hjust = 0.5),
                axis.title = element_text(size = 20,color ="black"),
                axis.text = element_text(size=20,color = "black"),
                panel.grid.minor.y = element_blank(),
                panel.grid.minor.x = element_blank(),
                legend.text = element_text(size=8),
                legend.title=element_blank(),
                axis.text.x = element_text(angle = 45, hjust=1, vjust=1)
)

count_table <- table(Immune@meta.data[,'group'], Immune@meta.data[,'mbcelltype'])
count_mtx <- as.data.frame.matrix(count_table)
count_mtx$cluster <- rownames(count_mtx)
melt_mtx <- melt(count_mtx)
melt_mtx$cluster <- as.factor(melt_mtx$cluster)

cluster_size <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)

if("0" %in% cluster_size$cluster){
  sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
}else{
  sorted_labels <- paste(cluster_size$cluster[order(cluster_size$value)])
}

cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)

colnames(melt_mtx)[2] <- "dataset"
p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") +
  theme_bw() + xlab("Cells per cluster") + ylab("") + mytheme

p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
  geom_bar(position="fill", stat="identity",) + theme_bw() +
  scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.9))+
  ylab(paste0("Fraction of cells in each ",tolower('cluster'))) + xlab('cluster') +
  theme(legend.position="top") + guides(fill = guide_legend(title ='cluster')) +
  scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme
wrap_plots(ncol = 2,p2,p1,widths =10)

melt_mtx$cluster <- factor(melt_mtx$cluster,levels=c("adjnromal","primary","metastasis"))

pdf('20230508_Immune_fraction.pdf',width = 5,height=7)
p2+scale_color_npg()
dev.off()
