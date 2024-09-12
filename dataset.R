rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(harmony)
##癌旁3##
dir_name <- c('GSM4710706','GSM4710707','GSM4710708')
adj <- list()
for (i in 1:length(dir_name)){
  dir.10x <- paste0("GSE155698_RAW/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  adj[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], 
                             min.cells = 3, min.features = 200)
}
print(i)
for (i in 1:length(adj)){
  sce <- adj[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce,pattern = "^MT-")
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce,pattern = "^RP[SL]")
  adj[[i]] <- sce
  rm(sce)
}
VlnPlot(adj[[2]],
        features = c("nFeature_RNA","nCount_RNA",
                     "percent.mt","percent.Ribo"),
        ncol = 4)
adata <- merge(adj[[1]],y = c(adj[[2]],adj[[3]]))
head(adata@meta.data)
adata@meta.data[["group"]] <- "Adjnorm"
#原位瘤的数据11+15#
dir_name <- c('GSM4679532','GSM4679533','GSM4679534','GSM4679535',
              'GSM4679535','GSM4679536','GSM4679537','GSM4679538',
              'GSM4679539','GSM4679540','GSM4679541')
p1 <- list()
for (i in 1:length(dir_name)){
  dir.10x <- paste0("GSE154778_RAW/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  p1[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], 
                             min.cells = 3, min.features = 200)
}
print(i)
dir_name <- c('GSM4710689','GSM4710690','GSM4710691','GSM4710692',
              'GSM4710693','GSM4710694','GSM4710695','GSM4710696',
              'GSM4710697','GSM4710698','GSM4710699','GSM4710700',
              'GSM4710701','GSM4710702','GSM4710704')
p2 <- list()
for (i in 1:length(dir_name)){
  dir.10x <- paste0("GSE155698_RAW/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  p2[[i]]=CreateSeuratObject(counts = my.data, project = dir_name[i], 
                             min.cells = 3, min.features = 200)
}
print(i)
plist <- list()
plist <- c(p1,p2)
for (i in 1:length(plist)){
  sce <- plist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce,pattern = "^MT-")
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce,pattern = "^RP[SL]")
  plist[[i]] <- sce
  rm(sce)
}
VlnPlot(plist[[5]],
        features = c("nFeature_RNA","nCount_RNA",
                     "percent.mt","percent.Ribo"),
        ncol = 4)
pdata <- merge(plist[[1]],y = c(plist[[2]], plist[[3]], plist[[4]],
                                plist[[5]], plist[[6]], plist[[7]],
                                plist[[8]], plist[[9]], plist[[10]],
                                plist[[11]],plist[[12]],plist[[13]],
                                plist[[14]],plist[[15]],plist[[16]],
                                plist[[17]],plist[[18]], plist[[19]], 
                                plist[[20]], plist[[21]], plist[[22]], 
                                plist[[23]], plist[[24]], plist[[25]], 
                                plist[[26]]))
head(pdata@meta.data)
pdata@meta.data[["group"]] <- "Primary"

#转移瘤的数据1+3+6#
a <- read.csv("PM_umiCounts_aboveBackground.csv")
row.names(a) <- a$CellId
a <- a[,-1]
b <- t(a)
m1 <- CreateSeuratObject(counts = b, project = "GSM4730268", 
                         min.cells = 3, min.features = 200)
dir_name <- c('GSM4730265','GSM4730266','GSM4730267')
m2 <- list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("GSE156405_RAW/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  m2[[i]]=CreateSeuratObject(counts = my.data, 
                             project = dir_name[i],
                             min.cells = 3, min.features = 200)
}

dir_name <- c('GSM4679542','GSM4679543','GSM4679544','GSM4679545',
              'GSM4679546','GSM4679547')
m3 <- list()
for (i in 1:length(dir_name)){
  dir.10x <-  paste0("GSE154778_RAW/",dir_name[i])
  my.data <- Read10X(data.dir = dir.10x) 
  m3[[i]]=CreateSeuratObject(counts = my.data, 
                             project = dir_name[i], 
                             min.cells = 3, min.features = 200)
}
print(i)
mlist <- list()
mlist <- c(m1,m2,m3)
for (i in 1:length(mlist)){
  sce <- mlist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce,pattern = "^MT-")
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce,pattern = "^RP[SL]")
  mlist[[i]] <- sce
  rm(sce)
}
VlnPlot(mlist[[5]],
        features = c("nFeature_RNA","nCount_RNA",
                     "percent.mt","percent.Ribo"),
        ncol = 4)
mdata <- merge(mlist[[1]],y = c(mlist[[2]], mlist[[3]], mlist[[4]],
                                mlist[[5]], mlist[[6]], mlist[[7]],
                                mlist[[8]], mlist[[9]], mlist[[10]]))
mdata@meta.data[["group"]] <- "Metastasis"
head(mdata@meta.data)
readdata <- merge(pdata,mdata,adata)
save(readdata, file = "readdata.RData")





##数据整合harmony##
rm(list = ls())
load("readdata.RData")
if(!require(harmony))devtools::install_github("immunogenomics/harmony")
library(harmony)
hdata <- readdata
hdata <-  hdata%>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData()
view(hdata)
hdata <- RunPCA(hdata, npcs = 50, verbose = F)
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = hdata, reduction = "pca", pt.size = .1, 
              group.by = "orig.ident")
p2 <- DimPlot(object = hdata, reduction = "pca", pt.size = .1, 
              group.by = "group")
p1+p2
p3 <- VlnPlot(object = hdata, features = "PC_1", pt.size = .1, 
              group.by = "orig.ident")
p4 <- VlnPlot(object = hdata, features = "PC_1", pt.size = .1, 
              group.by = "group")
p3|p4

#####run 到PCA再进行harmony，相当于降维########
hdata <- hdata %>% 
  RunHarmony("orig.ident",plot_convergence = TRUE,
             max.iter.harmony = 20)
hdata@reductions[["harmony"]]
harmony_embeddings <- Embeddings(hdata, 'harmony')
view(harmony_embeddings)
options(repr.plot.height = 5, repr.plot.width = 12)
p5 <- DimPlot(object = hdata, reduction = "harmony", pt.size = .1, 
              group.by = "orig.ident")
p6 <- VlnPlot(object = hdata, features = "harmony_1", pt.size = .1, 
              group.by = "orig.ident")
p7 <- DimPlot(object = hdata, reduction = "harmony", pt.size = .1, 
              group.by = "group")
p8 <- VlnPlot(object = hdata, features = "harmony_1", pt.size = .1, 
              group.by = "group")

hdata <- hdata %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

hdata <- hdata %>% 
  RunTSNE(reduction = "harmony", dims = 1:30)
TSNEPlot(object = hdata, pt.size = 0.5) 
TSNEPlot(object = hdata, pt.size = 0.5, label = TRUE,split.by='group') 

p5 <- DimPlot(hdata, reduction = "tsne", group.by = "orig.ident", pt.size=0.5)+
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),axis.text = element_blank()
  )
p5
save(hdata, file = "hdata.RData")


###降维###
ncol(hdata)
head(hdata@meta.data,5)
VlnPlot(readdata,
        features = c("nFeature_RNA",
                     "nCount_RNA",
                     "percent.mt","percent.Ribo"),
        ncol = 4,
        pt.size = 0.1) 
Idents(hdata)
Idents(hdata) <- "group"
plot1 <- FeatureScatter(hdata, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(hdata, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2
ggsave("nCount.pdf", p, width = 18, height = 12)


fdata <- subset(hdata,subset = nFeature_RNA > 200 & 
                  nFeature_RNA < 8000 & percent.mt < 25)
ncol(as.data.frame(fdata@assays[["RNA"]]@counts))
hdata <- fdata
hdata <- FindVariableFeatures(hdata, 
                              selection.method = "vst", 
                              nfeatures = 2000)
top10 <- head(VariableFeatures(hdata), 10)
top10
top20 <- head(VariableFeatures(hdata), 20)
top20
plot3 <- VariableFeaturePlot(object = hdata)
plot4 <- LabelPoints(plot = plot3, 
                     points = top10, 
                     repel = T)
plot3 + plot4
p <- plot3 + plot4
ggsave("volcano.pdf", p, width = 18, height = 12)

hdata <- ScaleData(hdata, features = rownames(hdata))
head(hdata[["RNA"]]@scale.data)
hdata <- RunPCA(hdata, 
                features = VariableFeatures(object = hdata))
print(hdata[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(hdata,dims = 1:4, reduction = "pca")
DimPlot(hdata, reduction = "pca")
DimHeatmap(hdata, dims = 1:30, cells = 500, balanced = TRUE)

hdata <- JackStraw(hdata,num.replicate = 100)
hdata <- ScoreJackStraw(hdata, dims = 1:20)
JackStrawPlot(hdata,dims = 1:20)
ElbowPlot(hdata,ndims = 30)

hdata <- FindNeighbors(hdata, dims = 1:26)
hdata <- FindClusters(hdata, resolution = 0.5)

hdata <- RunUMAP(hdata, dims = 1:26)
hdata <- RunTSNE(hdata,dims = 1:26,check_duplicates = FALSE)
DimPlot(hdata,reduction = "umap",label = T)
DimPlot(hdata, reduction = "tsne",label = T)
mdata <- hdata
save(mdata, file = "mdata.RData")


###找marker，注释###
markers <- FindAllMarkers(hdata, 
                          only.pos = TRUE, 
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

top10 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10,wt = avg_log2FC)
DoHeatmap(hdata, features = top10$gene) + NoLegend()
#DotPlot(hdata,features = top10$gene) + coord_flip()
save(top10,file = 'top10.RData')
save(markers,file = 'markers.RData')


cellmarker=c(
  "EPCAM",#上皮细胞 epithelial
  "PECAM1",#内皮细胞 endothelial
  "COL3A1",#成纤维细胞 fibroblasts
  "CD163","AIF1",#髓系细胞 myeloid
  "CD79A",#B细胞
  "JCHAIN",#浆细胞 plasma) cell
  "CD3D","CD8A","CD4",#T细胞
  "GNLY","NKG7",#NK细胞
  "PTPRC"#免疫细胞
)
VlnPlot(hdata,features = cellmarker,pt.size = 0,stack = T)

##singleR##
if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager") 
BiocManager::install(c("SingleR","celldex"),ask = F,update = F)
library(SingleR)
library(celldex)
library(pheatmap)
load("hpca.se.RData")
meta <- hdata@meta.data
head(meta)
plot1 <- DimPlot(hdata, reduction = "umap", label = TRUE)
plot2<-DimPlot(hdata, reduction = "tsne",label = TRUE)
plot1 + plot2

data_SingleR <- GetAssayData(hdata, slot="data") 
clusters <- hdata@meta.data$seurat_clusters
data.hesc <- SingleR(test = data_SingleR, ref = hpca.se, 
                     labels = hpca.se$label.main) 
view(data.hesc)
save(data.hesc,file = 'data.hesc.RData')
table(data.hesc$labels,meta$seurat_clusters)
tpye <- as.data.frame(table(data.hesc$labels,meta$seurat_clusters))
save(type,file = 'celltype.RData')

hdata@meta.data$labels <-data.hesc$labels
print(DimPlot(hdata, group.by = c("seurat_clusters", "labels"),
              reduction = "tsne"))


load('top10.RData')
c22 <- top10[221:230,]
VlnPlot(hdata,features = c22$gene,stack = T,pt.size = 0)
c28 <- top10[281:290,]
VlnPlot(hdata,features = c28$gene,stack = T,pt.size = 0)
c34 <- top10[341:350,]
VlnPlot(hdata,features = c34$gene,stack = T,pt.size = 0)
c35 <- top10[351:360,]
VlnPlot(hdata,features = c35$gene,stack = T,pt.size = 0)
##不确定的两群C31和C36##
c31 <- top10[311:320,]
VlnPlot(hdata,features = c31$gene,stack = T,pt.size = 0)
c36 <- top10[361:370,]
VlnPlot(hdata,features = c36$gene,stack = T,pt.size = 0)

FeaturePlot(hdata,features = c("GNLY","NKG7"))

hdata@meta.data$celltype <- recode(hdata@meta.data$seurat_clusters,
                                   '0'='Immune_cell',
                                   '1'='Immune_cell',
                                   '2'="Immune_cell",
                                   '3'='Epithelial_cell',
                                   '4'='Immune_cell',
                                   '5'='Immune_cell',
                                   '6'='Fibroblasts',
                                   '7'='Epithelial_cell',
                                   '8'='Epithelial_cell',
                                   '9'='Immune_cell',
                                   '10'='Epithelial_cell',
                                   '11'='Epithelial_cell',
                                   '12'='Epithelial_cell',
                                   '13'='Epithelial_cell',
                                   '14'='Epithelial_cell',
                                   '15'='Immune_cell',
                                   '16'='Immune_cell',
                                   '17'='Fibroblasts',
                                   '18'='Epithelial_cell',
                                   '19'='Immune_cell',
                                   '20'='Epithelial_cell',
                                   '21'='Immune_cell',
                                   '22'='Immune_cell',
                                   '23'='Endothelial_cell',
                                   '24'='Epithelial_cell',
                                   '25'='Immune_cell',
                                   '26'='Epithelial_cell',
                                   '27'='Epithelial_cell',
                                   '28'='Immune_cell',
                                   '29'='Immune_cell',
                                   '30'='Fibroblasts',
                                   '31'='Epithelial_cell',
                                   '32'='Immune_cell',
                                   '33'='Epithelial_cell',
                                   '34'='Immune_cell',
                                   '35'='Immune_cell',
                                   '36'='Immune_cell',
                                   '37'='Immune_cell')
unique(hdata@meta.data$celltype)
save(hdata,file = 'zhushidata.RData')

DimPlot(hdata, group.by="celltype", reduction='tsne')
DimPlot(hdata, group.by="celltype", label=T, label.size=6, 
        reduction='umap')
##ggplot2可视化##
phe <- data.frame(cell=rownames(hdata@meta.data),
               cluster =clusters,
               celltype=hdata@meta.data$celltype)
tsne_pos <- Embeddings(hdata,'tsne')
dat <- cbind(tsne_pos,phe)
head(dat)
p <- ggplot(dat,aes(x=tSNE_1,y=tSNE_2,color=celltype))+
  geom_point(size=0.95)
p
p <- p+stat_ellipse(data=dat,aes(x=tSNE_1,y=tSNE_2,
                              fill=celltype,color=celltype),
                 geom = "polygon",alpha=0.2,level=0.9,
                 type="t",linetype = 2,show.legend = F)+
  coord_fixed()
p
theme= theme(panel.grid =element_blank()) +   ## 删去网格
  theme(panel.border = element_blank(),panel.background = element_blank()) +   ## 删去外层边框
  theme(axis.line = element_line(size=1, colour = "black"))
p <- p+theme+guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("tsne_celltype.pdf", p, width = 18, height = 12)







##寻找差异基因##
celldata <- hdata
table(celldata@meta.data$celltype)
names(celldata@meta.data)
unique(celldata$group)
DimPlot(celldata)
DimPlot(celldata,split.by = 'group',label = T)
#idents重要性##
Idents(celldata)
celldata$celltype.group <- paste(celldata$celltype, 
                                 celldata$group, sep = "_")
unique(celldata$celltype.group)
celldata$celltype <- Idents(celldata)
Idents(celldata) <- "celltype.group"
##单个细胞##
##adjnorm和primary##
f <- FindMarkers(celldata,
                 ident.1 = 'Fibroblasts_Adjnorm',
                 ident.2 = 'Fibroblasts_Primary',
                 verbose = T, test.use = 'DESeq2',min.pct = 0.1,slot='conuts')
head(f)
write.csv(f,"AP_deg.Fibroblasts.CSV")
i <- FindMarkers(celldata,
                 ident.1 = 'Immune_cell_Adjnorm',
                 ident.2 = 'Immune_cell_Primary',
                 verbose = T, test.use = 'wilcox',min.pct = 0.1)
head(i)
write.csv(i,"AP_deg.Immune_cell.CSV")
epi <- FindMarkers(celldata,
                 ident.1 = 'Epithelial_cell_Adjnorm',
                 ident.2 = 'Epithelial_cell_Primary',
                 verbose = T, test.use = 'wilcox',min.pct = 0.1)
head(epi)
write.csv(epi,"AP_deg.Epithelial_cell.CSV")
endo <- FindMarkers(celldata,
                   ident.1 = 'Endothelial_cell_Adjnorm',
                   ident.2 = 'Endothelial_cell_Primary',
                   verbose = T, test.use = 'wilcox',min.pct = 0.1)
head(endo)
write.csv(endo,"AP_deg.Endothelial_cell.CSV")
##primary和metastasis##
ff <- FindMarkers(celldata,
                  ident.1 = 'Fibroblasts_Primary',
                  ident.2 = 'Fibroblasts_Metastasis',
                  verbose = T, test.use = 'DESeq2',min.pct = 0.1)
head(ff)
write.csv(ff,"PM_deg.fibroblasts.CSV")
ii <- FindMarkers(celldata,
                 ident.1 = 'Immune_cell_Primary',
                 ident.2 = 'Immune_cell_Metastasis',
                 verbose = T, test.use = 'wilcox',min.pct = 0.1)
head(ii)
write.csv(ii,"PM_deg.Immune_cell.CSV")
ee <- FindMarkers(celldata,
                   ident.1 = 'Epithelial_cell_Primary',
                   ident.2 = 'Epithelial_cell_Metastasis',
                   verbose = T, test.use = 'wilcox',min.pct = 0.1)
head(ee)
write.csv(ee,"PM_deg.Epithelial_cell.CSV")
eded <- FindMarkers(celldata,
                    ident.1 = 'Endothelial_cell_Primary',
                    ident.2 = 'Endothelial_cell_Metastasis',
                    verbose = T, test.use = 'wilcox',min.pct = 0.1)
head(eded)
write.csv(eded,"PM_deg.Endothelial_cell.CSV")

##两组间循环计算##
cellfordeg<-levels(celldata$celltype)
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(celldata, ident.1 = paste0(cellfordeg[i],"_Adjnorm"), 
                         ident.2 = paste0(cellfordeg[i],"_Primary"), 
                         verbose = T)
  write.csv(CELLDEG,paste0(cellfordeg[i],".CSV"))
}

##EnhancedVolcano##
##BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
deg.Fibroblast <- read.csv("PM_deg.Fibroblasts.CSV",row.names = 1)
deg.Immune_cell <- read.csv("PM_deg.Immune_cell.CSV",row.names = 1)
deg.Endothelial_cell <- read.csv("PM_deg.Endothelial_cell.CSV",row.names = 1)
deg.Epithelial_cell <- read.csv("PM_deg.Epithelial_cell.CSV",row.names = 1)
EnhancedVolcano(deg.Fibroblast,
                lab = rownames(deg.Fibroblast),
                x = 'avg_log2FC',y = 'p_val_adj',
                pCutoff = 0.05,FCcutoff = 0.5,
                pointSize = 3,labSize = 5,
                title = 'deg.Fibroblast')

deg.Fibroblast$celltype <- "Fibroblast"
deg.Endothelial_cell$celltype <- "Endothelial_cell"
deg.Epithelial_cell$celltype <- "Epithelial_cell"
deg.Immune_cell$celltype <- "Immune_cell"
diff.cell <- rbind(deg.Immune_cell,deg.Epithelial_cell,
                   deg.Fibroblast,deg.Endothelial_cell)
diff.cell$gene <- rownames(diff.cell)
head(diff.cell)

diff.cell$change <- ifelse(diff.cell$p_val_adj<0.01,
                           "adjust P-val<0.01",
                           "adjust p-val>=0.01")
top10_Fibroblast <- filter(diff.cell,celltype=="Fibroblast") %>%
  distinct(gene,.keep_all=T)%>%
  top_n(n = 10,wt = abs(avg_log2FC))
top10_Endothelial_cell <- filter(diff.cell,celltype=="Endothelial_cell") %>%
  distinct(gene,.keep_all=T)%>%
  top_n(n = 10,wt = abs(avg_log2FC))
top10_Epithelial_cell <- filter(diff.cell,celltype=="Epithelial_cell") %>%
  distinct(gene,.keep_all=T)%>%
  top_n(n = 10,wt = abs(avg_log2FC))
top10_Immune_cell <- filter(diff.cell,celltype=="Immune_cell") %>%
  distinct(gene,.keep_all=T)%>%
  top_n(n = 10,wt = abs(avg_log2FC))
top10 <-rbind(top10_Epithelial_cell,
              top10_Immune_cell,
              top10_Fibroblast,
              top10_Endothelial_cell)
write.csv(top10,"PM_deg.top10.CSV")
diff.cell$size <-case_when(!(diff.cell$gene %in% 
                               top10$gene)~ 1,
                           diff.cell$gene %in% top10$gene~2)
##提取不显示标题的数据，将其作图
dt<-filter(diff.cell,size==1)
##绘制背景柱和散点，因为横轴是因子，所以后面的作图中要注意
dfbar<-data.frame(x=c("Immune_cell","Epithelial_cell",
                      "Fibroblast","Endothelial_cell"),
                  y=c(5,5,5,5))
dfbar1<-data.frame(x=c("Immune_cell","Epithelial_cell",
                       "Fibroblast","Endothelial_cell"),
                   y=c(-5.5,-5.5,-5.5,-5.5))
p<- ggplot()+
  geom_col(data = dfbar,mapping=aes(x=x,y=y),
           fill = "#dcdcdc",alpha=0.6)+
  geom_col(data=dfbar1,mapping=aes(x=x,y=y),
           fill = "#dcdcdc",alpha =0.6)+
  geom_jitter(data = dt,
              aes(x=celltype,y= avg_log2FC,color =change),
              size=0.85,width=0.4)+
  geom_jitter(data = top10,
              aes(x= celltype,y= avg_log2FC, color = change),
              size=1.5,width=0.4)
p
##添加显著基因名##
dfcol<-data.frame(x=c("Immune_cell","Epithelial_cell",
                      "Fibroblast","Endothelial_cell"),
                  y=0,
                  label=c("Immune_cell","Epithelial_cell",
                          "Fibroblast","Endothelial_cell"))
mycol <- c("#DAA520","#3CB371","#0000CD","#FF44AA")
library(ggrepel)
p1 <- p + geom_tile(data = dfcol,aes(x=x,y=y),
                    height=0.4,
                    color = "black",
                    fill= mycol,
                    alpha=0.5,
                    show.legend=F)+
  geom_text_repel(
    data = top10,
    aes(x = celltype,y = avg_log2FC,label = gene),
    size=3,
    max.overlaps = 20,
  )
p1
p2<-p1+
  scale_color_manual(name=NULL,
                     values = c("red","grey"))+
  labs(x="celltype",y="avg_log2FC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size=5,
            color ="white")+
  theme_minimal()+
  theme(
    axis.title = element_text(size =12,
                              color ="black",
                              face="bold"),
    axis.line.y = element_line(color = "black",
                               size=0.8),
    axis.line.x=element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification=c(1,0),
    legend.text=element_text(size=12)
  )
p2
ggsave("EnhancedVolcano_PM.pdf", p2, width = 10, height = 8)
##常规火山图##








##富集分析##
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db) ##加载人类
library(tidyverse)
library(gplots)
library(ggplot2)
deg.immune <- read.csv("PM_deg.Immune_cell.CSV",row.names = 1)
sig_dge.all <- subset(deg.immune, p_val_adj<0.05&
                        abs(avg_log2FC)>0.15)
View(sig_dge.all)
ego_ALL <- enrichGO(gene          = row.names(sig_dge.all),
                    #universe     = row.names(dge.celltype),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",  #设置为ALL时BP, CC, MF都计算
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)
View(ego_all)
ego_CC <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_cc <- data.frame(ego_CC)
ego_MF <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_mf <- data.frame(ego_MF)
ego_BP <- enrichGO(gene          = row.names(sig_dge.all),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05) 
ego_bp <- data.frame(ego_BP)
p_BP <- barplot(ego_BP,showCategory = 10) + 
  ggtitle("barplot for Biological process")
p_CC <- barplot(ego_CC,showCategory = 10) + 
  ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10) + 
  ggtitle("barplot for Molecular function")
plotc <- p_BP/p_CC/p_MF
plotc
##ggplot2##
ego_bp <- ego_bp[order(ego_bp$p.adjust),]
ego_bp_top30 <- ego_bp[1 : 30,]
ggplot(data=ego_bp_top30, aes(x=Description,y=Count)) + 
  geom_bar(stat="identity", width=0.8,fill='salmon1') + 
  coord_flip() +  xlab("GO term") + ylab("Num of Genes") + 
  theme_bw()








##KEGG##
genelist <- bitr(row.names(sig_dge.all), fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
p1 <- barplot(ekegg, showCategory=20)
p2 <- dotplot(ekegg, showCategory=20)
plotc = p1/p2
plotc
##GSVA##
top10gene <- deg.immune %>% 
  top_n(n = 10,wt = avg_log2FC)%>%rownames()
buttom10gene <- deg.immune %>% 
  top_n(n = -10,wt = avg_log2FC)%>%rownames()
geneset <- as.data.frame(cbind(top10gene,buttom10gene))
load("hdata.RData")
exp <- as.matrix(hdata@assays[["RNA"]]@data)
dim(exp)
if(!require("GSVA")) install.packages("GSVA",update = F,ask = F)
library(GSVA)
gsva.result <- GSVA::gsva(exp,geneset,kcdf = "Gaussian")
mymatrix <- t(gsva.result)
mymatrix <- cbind(mymatrix,as.data.frame(hdata$labels))
colnames(mymatrix)[ncol(mymatrix)] <- "celltype"
head(mymatrix)
