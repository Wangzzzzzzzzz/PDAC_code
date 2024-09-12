if (!require("harmony", quietly = TRUE))
     install.packages("harmony")
install.packages(c("SeuratObject", "Seurat"))
library(SeuratObject)
library(Seurat)

library(harmony)
library(sctransform) 
library(future)

scRNA<- merge(obj.list[[1]], y=obj.list[-1])
nbrOfWorkers()
plan(multisession, workers=2)

dim(scRNA)
table(scRNA@meta.data$orig.ident)

scRNA_harmony <- NormalizeData(scRNA,verbose = T) 
scRNA_harmony <- FindVariableFeatures(scRNA_harmony,verbose = T)  
scRNA_harmony <-  ScaleData(scRNA_harmony,verbose = T)
scRNA_harmony <-  RunPCA(scRNA_harmony)

system.time({scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")})
obj.combined <- scRNA_harmony
ElbowPlot(obj.combined, ndims = 50)
saveRDS(obj.combined,'20230427_large_harmonied_LogNorm.rds')

setwd("/home/shpcv2_kvcka/data/20230426_PDAC")
getwd()
obj.combined <- readRDS('20230427_large_harmonied_LogNorm.rds')
obj.combined <- RunUMAP(obj.combined, reduction = "harmony", dims = 1:40)
DimPlot(obj.combined)

obj.combined <- FindNeighbors(obj.combined, reduction = "harmony", dims = 1:40)

obj.combined <- FindClusters(obj.combined,resolution = 0.8)
DimPlot(obj.combined)

library(cowplot)
library(ggplot2)

sel.clust = "RNA_snn_res.0.8"
obj.combined <- SetIdent(obj.combined, value = sel.clust)
table(obj.combined@active.ident)

pdf('20230427_obj.combined_umap.pdf',width = 11,height = 9)
DimPlot(obj.combined, reduction = "umap", group.by = "GSM",pt.size = 0.01)
DimPlot(obj.combined, reduction = "umap", label = TRUE,pt.size = 0.01)
dev.off()


FeaturePlot(obj.combined,features = c("PCNA","MKI67"),reduction = "umap")
RidgePlot(obj.combined, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

cc_gene=unlist(cc.genes)
s.genes=cc.genes$s.genes

library(Hmisc)

g2m.genes=cc.genes$g2m.genes

obj.combined <- CellCycleScoring(obj.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj.combined@meta.data  %>% ggplot(aes(x=S.Score,y=G2M.Score))+geom_point(aes(color=Phase))+theme_minimal()

pdf('20230427_umap_groupby_phases.pdf')
DimPlot(obj.combined,group.by = 'Phase')
dev.off()
table(obj.combined$Phase)
saveRDS(obj.combined,'20230427_umaped_obj.combined.rds')
