setwd("/home/shpcv2_kvcka/data/20230426_PDAC/Annotation")
getwd()

sel.clust = "RNA_snn_res.0.8"
obj.combined <- SetIdent(obj.combined, value = sel.clust)
、
pdf('20230428_bigcelltype_featureplot.pdf',height =7,width = 10)
FeaturePlot(obj.combined,features = 'PTPRC')    
FeaturePlot(obj.combined,features = c('KRT18','EPCAM')) 
FeaturePlot(obj.combined,features = c('PECAM1','COL3A1'))  
dev.off()


sel.clust = "RNA_snn_res.0.8"
obj.combined <- SetIdent(obj.combined, value = sel.clust)
table(obj.combined@active.ident)

mysub <- FindSubCluster(
  obj.combined,
  graph.name='RNA_snn', 
  cluster='9', 
  subcluster.name = "sub.cluster",
  resolution = 0.5,
  algorithm = 1                           
)
Idents(mysub) <- 'sub.cluster'

pdf(file="2023428_obj.combined_sub.cluster.pdf",width=15,height = 15)
DimPlot(mysub,label = T,label.size = 2)
dev.off()

obj.combined <- mysub
rm(mysub)
gc()

pdf('20230427_umap_markers_dotplot.pdf',width = 10,height =8)
DotPlot(obj.combined,assay = "RNA",cols = rep(c("darkblue", "darkred",'darkorange','black','darkgreen'),8),features = genes_to_check, split.by = "RNA_snn_res.0.8") + coord_flip()+ RotatedAxis()
dev.off()
pdf('20230428_umap_sub.clusters.pdf',width = 10,height =8)
DimPlot(obj.combined,label = T,label.size =2)
dev.off()

obj.combined <- SetIdent(obj.combined, value = "sub.cluster")
table(obj.combined@active.ident)
new.cluster.ids <- c( "5","0","17","9_0","1",
                      "9_2","7","11", "15","9_1",
                      "2","22","18","24","29","9_4",
                      "19","12","27","31","4","20",
                      "3","8","6","16","25","9_3",
                      "21","10","14","28","13","9_5",
                      "32","9_8","9_9","34","33","36",
                      "30","23","26","35","9_7","9_6" )
new.cluster.ids <- c( "Epithelial_cells","Epithelial_cells","Epithelial_cells","Epithelial_cells","Immune_cells",
                      "Epithelial_cells","Immune_cells","Epithelial_cells", "Immune_cells","Epithelial_cells",
                      "Immune_cells","Epithelial_cells","Immune_cells","Epithelial_cells","Stromal_cells","Immune_cells",
                      "Stromal_cells","Immune_cells","Immune_cells","Epithelial_cells","Immune_cells","Epithelial_cells",
                      "Immune_cells","Immune_cells","Stromal_cells","Stromal_cells","Immune_cells","Immune_cells",
                      "Epithelial_cells","Stromal_cells","Immune_cells","Epithelial_cells","Immune_cells","Epithelial_cells",
                      "Epithelial_cells","Epithelial_cells","Epithelial_cells","Stromal_cells","Immune_cells","Immune_cells",
                      "Immune_cells","Immune_cells","Epithelial_cells","Immune_cells","Immune_cells","Epithelial_cells" )
names(new.cluster.ids) <- levels(obj.combined)
obj.combined <- RenameIdents(obj.combined, new.cluster.ids)
table(obj.combined@active.ident)
pdf('20230428_umap_annotations.pdf',width = 9,height = 7)
DimPlot(obj.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()
obj.combined$bigcelltype <- obj.combined@active.ident
saveRDS(obj.combined,'20230428_obj.combined_bigcelltype_annotated.rds')

Epithelial <- subset(obj.combined,idents = 'Epithelial_cells')
Immune <- subset(obj.combined,idents = 'Immune_cells')
Stromal<- subset(obj.combined,idents = 'Stromal_cells')


library(harmony)
Epithelial@active.assay
Epithelial_harmony <- SCTransform(object = Epithelial,vars.to.regress = "percent.MT")

system.time({Epithelial_harmony <- RunHarmony(Epithelial_harmony, group.by.vars = "orig.ident")})
Epithelial <- Epithelial_harmony
rm(Epithelial_harmony)
DimPlot(Epithelial,reduction = 'harmony',group.by = 'orig.ident')
ElbowPlot(Epithelial, ndims = 50)
saveRDS(Epithelial,'/home/shpcv2_kvcka/data/20230426_PDAC/Epithelial/Epithelial_harmonied.rds')

library(harmony)
Immune@active.assay
Immune_harmony <- SCTransform(object = Immune,vars.to.regress = "percent.MT")
system.time({Immune_harmony <- RunHarmony(Immune_harmony, group.by.vars = "orig.ident")})
Immune <- Immune_harmony
rm(Immune_harmony)
DimPlot(Immune,reduction = 'harmony',group.by = 'orig.ident')
ElbowPlot(Immune, ndims = 50)
saveRDS(Immune,'/home/shpcv2_kvcka/data/20230426_PDAC/Immune/Immune_harmonied.rds')


library(harmony)
Stromal@active.assay
Stromal_harmony <- SCTransform(object = Stromal,vars.to.regress = "percent.MT")
system.time({Stromal_harmony <- RunHarmony(Stromal_harmony, group.by.vars = "orig.ident")})
Stromal <- Stromal_harmony
rm(Stromal_harmony)
DimPlot(Stromal,reduction = 'harmony',group.by = 'orig.ident')
ElbowPlot(Stromal, ndims = 50)
saveRDS(Stromal,'/home/shpcv2_kvcka/data/20230426_PDAC/Stromal/Stromal_harmonied.rds')




obj.combined <- readRDS('20230428_obj.combined_bigcelltype_annotated.rds')



DefaultAssay(Immune) <- "SCT"
VlnPlot(Immune,features = 'PTPRC',group.by = 'sub.cluster')
 
pdf('20230427_umap_markers_dotplot.pdf',width = 10,height =8)
DotPlot(obj.combined,assay = "RNA",cols = rep(c("darkblue", "darkred",'darkorange','black','darkgreen'),8),features = genes_to_check, split.by = "RNA_snn_res.0.8") + coord_flip()+ RotatedAxis()
dev.off()

color <- c('lightgrey','darkred')
pdf('20230507_SLC12A5_FeaturePlot.pdf')
FeaturePlot(obj.combined,order=T,reduction='umap',features = c('SLC12A5'),ncol = 1,cols=color,pt.size=0.5,label.size = 20)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+
  theme(plot.title = element_text(size = 45))
dev.off()
pdf('20230430_KRT18_FeaturePlot.pdf')
FeaturePlot(obj.combined,order=F,reduction='umap',features = c('KRT18'),ncol = 1,cols=color,pt.size=0.1,label.size = 20)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+
  theme(plot.title = element_text(size = 45))
dev.off()
pdf('20230430_EPCAM_FeaturePlot.pdf')
FeaturePlot(obj.combined,order=F,reduction='umap',features = c('EPCAM'),ncol = 1,cols=color,pt.size=0.1,label.size = 20)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+
  theme(plot.title = element_text(size = 45))
dev.off()
pdf('20230430_PECAM1_FeaturePlot.pdf')
FeaturePlot(obj.combined,order=F,reduction='umap',features = c('PECAM1'),ncol = 1,cols=color,pt.size=0.1,label.size = 20)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+
  theme(plot.title = element_text(size = 45))
dev.off()
pdf('20230430_COL3A1_FeaturePlot.pdf')
FeaturePlot(obj.combined,order=F,reduction='umap',features = c('COL3A1'),ncol = 1,cols=color,pt.size=0.1,label.size = 20)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+
  theme(plot.title = element_text(size = 45))
dev.off()
pdf('20230430_PTPRC_FeaturePlot.pdf')
FeaturePlot(obj.combined,order=F,reduction='umap',features = c('PTPRC'),ncol = 1,cols=color,pt.size=0.1,label.size = 20)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+
  theme(plot.title = element_text(size = 45))
dev.off()
pdf('20230430_MKI67_FeaturePlot.pdf')
FeaturePlot(obj.combined,order=F,reduction='umap',features = c('MKI67'),ncol = 1,cols=color,pt.size=0.1,label.size = 20)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+
  theme(plot.title = element_text(size = 45))
dev.off()




obj.combined <- readRDS('20230428_obj.combined_bigcelltype_annotated.rds')
Idents(obj.combined) <- 'RNA_snn_res.0.8'
DimPlot(obj.combined,label = T)
macro <- subset(obj.combined,idents=c('1','33','7'))
macro_stomal <- subset(obj.combined,idents=c('1','33','7','6','16','10','19'))
saveRDS(macro,'/home/shpcv2_kvcka/data/20230426_PDAC/Macro/macro_raw.rds')
saveRDS(macro_stomal,'/home/shpcv2_kvcka/data/20230426_PDAC/Macro_Stromal/Macro_Stromal_raw.rds')




color <- c('lightgrey','darkred')
genes <- c('MAFB','ENPP2','NDFIP1','PLD3','ALDH1A1','TSPAN33','ANKH','SLC7A8','CD163','SLC12A5')
for (i in genes){
  print(i)
  pdf(paste0('/home/shpcv2_kvcka/data/20230426_PDAC/TCGA/20230504_',i,'_FeaturePlot.pdf'))
  p <- FeaturePlot(obj.combined,order=F,reduction='umap',features = i,ncol = 1,cols=color,pt.size=0.3,label.size = 20)+
    theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+
    theme(plot.title = element_text(size = 45))
  print(p)
  dev.off()
  
}

