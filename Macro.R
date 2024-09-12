setwd('/home/shpcv2_kvcka/data/20230426_PDAC/Macro')
getwd()
macro <- readRDS('macro_raw.rds')

library(harmony)
macro@active.assay
macro_harmony <- SCTransform(object = macro,vars.to.regress = "percent.MT")

system.time({macro_harmony <- RunHarmony(macro_harmony, group.by.vars = "orig.ident")})
macro <- macro_harmony
rm(macro_harmony)
DimPlot(macro,reduction = 'harmony',group.by = 'orig.ident')
ElbowPlot(macro, ndims = 50)


macro@graphs$graphname
macro <- RunTSNE(macro, reduction = "harmony", dims = 1:40)
macro <- FindNeighbors(macro, reduction = "harmony", dims = 1:40)
macro <- FindClusters(macro,resolution = 0.8,graph.name = "RNA_snn")
pdf('20230501_macro_tsne.pdf')
DimPlot(macro,reduction = 'tsne',label=T)
dev.off()
FeaturePlot(macro,features = 'SLC12A5',order=T,reduction = 'tsne',pt.size = 1)
new.cluster.ids <- c('SLC12A5_high_M2',"SLC12A5_low_M2","SLC12A5_high_M2", "SLC12A5_low_M2","M1-like TAM","M1-like TAM",
                     'M1-like TAM',"SLC12A5_low_M2","M1-like TAM","SLC12A5_high_M2","SLC12A5_low_M2",
                     "SLC12A5_low_M2","SLC12A5_low_M2","13")
table(macro@active.ident)
names(new.cluster.ids) <- levels(macro)
macro <- RenameIdents(macro, new.cluster.ids)
table(macro@active.ident)
macro <- subset(macro,ident=c('SLC12A5_high_M2','SLC12A5_low_M2','M1-like TAM'))
macro$SLC12A5_macro <- macro@active.ident
macro$SLC12A5 <- NULL
pdf('20230501_macro_annotated.pdf')
DimPlot(macro,reduction = 'tsne',label = T)
dev.off()
saveRDS(macro,'20230501_macro_annotated.rds')
#########
Idents(macro) <- 'RNA_snn_res.0.8'
TSNEPlot(macro,label=T)
new.cluster.ids <- c('SLC12A5_high_M2',"SLC12A5_low_M2","SLC12A5_low_M2", "SLC12A5_low_M2","M1-like TAM","M1-like TAM",
                     'M1-like TAM',"SLC12A5_low_M2","M1-like TAM","SLC12A5_high_M2","SLC12A5_low_M2",
                     "SLC12A5_low_M2","SLC12A5_low_M2")
table(macro@active.ident)
names(new.cluster.ids) <- levels(macro)
macro <- RenameIdents(macro, new.cluster.ids)
table(macro@active.ident)
macro$SLC12A5_macro2 <- macro@active.ident


library(ggplot2)
color <- c('lightgrey','darkred')
pdf('20230501_SLC12A5_in_macro.pdf')
FeaturePlot(macro,features = 'SLC12A5',ncol = 1,cols=color,pt.size=1,reduction = 'tsne',order=T)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+
  theme(plot.title = element_text(size = 35))
dev.off()

pdf('20230501_macro_annotated.pdf',height = 10,width = 10)
DimPlot(macro,ncol = 1,pt.size=1,reduction = 'tsne',order=T,label = T,label.size = 10)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+
  theme(plot.title = element_text(size = 35))+scale_color_npg()+theme(legend.position = 'none')
dev.off()


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

count_table <- table(macro@meta.data[,'group'], macro@meta.data[,'SLC12A5_macro'])
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

pdf('20230501_macro_fraction.pdf',width = 5,height=7)
p2+guides(fill = guide_legend( ncol = 1, byrow = TRUE))+scale_color_npg()
dev.off()


rm(list=ls())
gc()

macro <- readRDS('20230501_macro_annotated.rds')
#DefaultAssay(macro) <- "SCT"
M2.diff <- FindMarkers(macro, ident.1 = "SLC12A5_high_M2", ident.2 = "SLC12A5_low_M2", logfc.threshold = 0,verbose = T)
write.csv(file="20230501_M2.diff.csv",M2.diff)
M2.all <- M2.diff
M2.diff <- subset(M2.diff,p_val_adj<0.05)

low<-floor(range(M2.diff$avg_log2FC)[1]) 
high<-ceiling(range(M2.diff$avg_log2FC)[2]) 
pdf('20230501_DEG_M2.pdf',width = 10,height=10)
print(EnhancedVolcano(M2.diff,
                      title = 'SLC12A5_high_M2 versus SLC12A5_low_M2',
                      lab = rownames(M2.diff),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 0.2,
                      xlim = c(low, high))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
dev.off()

## KEGG pathway analysis
R.utils::setOption("clusterProfiler.download.method",'auto')

rt <- M2.diff
rt$id <- rownames(rt)
rt <- subset(rt,subset=p_val_adj<0.05)
rt <- subset(rt,subset=abs(avg_log2FC)>0.20)
rt <- rt[,c(-1,-3,-4,-5)]
genes=as.vector(rt[,2])  
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)   
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
kk.up <- enrichKEGG(gene         = out$entrezID,
                    organism     = 'hsa',
                    pvalueCutoff = 0.05,
                    qvalueCutoff =0.05)
head(kk.up)[,1:6]


kk.up = setReadable(kk.up,
                    OrgDb = "org.Hs.eg.db",
                    keyType = "ENTREZID")


kegg_up_dt <- as.data.frame(kk.up)
up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,]

write.csv(file="20230510_M2.diff_KEGG_.csv",kegg_up_dt)


##
kegg_up_dt$p <- kegg_up_dt$p.adjust
for (i in 1:nrow(kegg_up_dt)){
  print(i)
  kegg_up_dt$Description[i] <- substr(as.character(kegg_up_dt$Description[i]),1,(nchar(as.character(kegg_up_dt$Description[i]))-23))
}
kegg_up_dt <- kegg_up_dt[c(2,4,7,10,17,20,21,23,24,26,27,29,30,31,32,33,36,37,38,40,44,45,47,49,50,51,52,54,60,63,67,69,73,74,75,76,77),]

pdf(file="20230510_M2.diff_KEGG.pdf",width=9,height =5)
ggplot(kegg_up_dt,aes(Count,Description))+
  geom_bar(aes(y=reorder(Description,Count),x=Count,fill=p)
           ,stat='identity')+
  scale_fill_gradient(low="#FFCC33",high="#CC6666")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))
dev.off()

kegg_down_dt$p <- kegg_down_dt$p.adjust

rt <- M2.diff
rt$id <- rownames(rt)
rt <- subset(rt,subset=p_val_adj<0.05)
rt <- rt[,c(-1,-3,-4,-5)]
genes=as.vector(rt[,2])  
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
geneList= rt$avg_log2FC 
names(geneList)= rt$entrezID
geneList=sort(geneList,decreasing = T)
head(geneList)
kk_gse <- gseKEGG(geneList     = geneList,
                  organism     = 'hsa',
                  nPerm        = 1000,
                  minGSSize    = 10,
                  pvalueCutoff = 0.1,
                  verbose      = T)
head(kk_gse)[,1:6]

write.csv(file="20230501_SLC12A5HighM2.vs.SLC12A5LowM2.diff_GSEA_KEGG_result.csv",data.frame(kk_gse))
for (i in 1:length(kk_gse$Description)){
  print(i)
  kk_gse$Description[i] <- substr(as.character(kk_gse$Description[i]),1,(nchar(as.character(kk_gse$Description[i]))-23))
}
pdf(file="20230501_SLC12A5HighM2.vs.SLC12A5LowM2.diff_GSEA_KEGG.pdf",height = 17,width=15)
clusterProfiler::dotplot(kk_gse,showCategory=50,label_format=80) #showCategory显示的通路
dev.off()


save(kk_gse,M2.all,M2.diff,file = '20230502GSEA.Rdata')


load('20230502GSEA.Rdata')

library(GseaVis)
pdf(file="GSEAResult/20230510_GSEA_hsa04658.pdf",width=4,height=3)
gseaNb(object = kk_gse,
       geneSetID = 'hsa04658',
       newGsea = F)
dev.off()

macro <- readRDS('20230501_macro_annotated.rds')

library(monocle)
macro$celltype1 <- Idents(macro)
data <- as(as.matrix(macro@assays$RNA@counts), 'sparseMatrix')

p_data <- macro@meta.data
p_data$celltype <- macro@active.ident 
pd <- new('AnnotatedDataFrame', data = p_data) 

f_Data <- data.frame(gene_short_name = row.names(macro), row.names = row.names(macro))
fd <- new('AnnotatedDataFrame', data = f_Data)
identical(colnames(data),rownames(pd))
identical(rownames(data),rownames(fd))
data <- data[rownames(fd),]
HSMM_all <- newCellDataSet(data,      
                           phenoData = pd,
                           featureData = fd)


ordering_genes <- macro[["SCT"]]@var.features
HSMM_all <- setOrderingFilter(HSMM_all, ordering_genes)
HSMM_all <- estimateSizeFactors(HSMM_all)  
HSMM_all <- estimateDispersions(HSMM_all) 
HSMM_all <- reduceDimension(HSMM_all,
                            norm_method="none",
                            reduction_method="DDRTree",  
                            max_components=3,
                            scaling=TRUE,
                            verbose=TRUE,
                            pseudo_expr=0)
HSMM_all <- orderCells(HSMM_all) 
saveRDS(HSMM_all,'20230501_HSMM_all_SCT.rds')


library(monocle)
macro$celltype1 <- Idents(macro)
macro_primary <- subset(macro,subset=group=='primary')
data <- as(as.matrix(macro_primary@assays$RNA@counts), 'sparseMatrix')

p_data <- macro_primary@meta.data
p_data$celltype <- macro_primary@active.ident  
pd <- new('AnnotatedDataFrame', data = p_data)

f_Data <- data.frame(gene_short_name = row.names(macro_primary), row.names = row.names(macro_primary))
fd <- new('AnnotatedDataFrame', data = f_Data)
identical(colnames(data),rownames(pd))
identical(rownames(data),rownames(fd))
data <- data[rownames(fd),]
HSMM_primary <- newCellDataSet(data,       
                           phenoData = pd,
                           featureData = fd)


ordering_genes <- macro_primary[["RNA"]]@var.features
HSMM_primary <- setOrderingFilter(HSMM_primary, ordering_genes)
HSMM_primary <- estimateSizeFactors(HSMM_primary)  
HSMM_primary <- estimateDispersions(HSMM_primary)  
HSMM_primary <- reduceDimension(HSMM_primary,
                            norm_method="none",
                            reduction_method="DDRTree",  
                            max_components=3,
                            scaling=TRUE,
                            verbose=TRUE,
                            pseudo_expr=0)
HSMM_primary <- orderCells(HSMM_primary) 
saveRDS(HSMM_primary,'20230501_HSMM_primary_SCT.rds')
saveRDS(HSMM_primary,'20230501_HSMM_primary_Log.rds')

pdf('20230430_trajectory_state.pdf')
plot_cell_trajectory(HSMM_primary)
dev.off()

colour <- c('#E64B35FF','#4DBBD5FF','#00A087FF','#3C54BBFF','#F39B7FFF',
            '#8491B4FF','#91D1C2FF','#DC0000FF','#7E6148FF')
pdf('10_trajectory_celltype.pdf')
plot_cell_trajectory(HSMM_primary, color_by = 'SLC12A5_macro',order=F)+scale_color_npg()
dev.off()
# plot_cell_trajectory(HSMM_all, color_by = "Pseudotime")+ scale_color_manual(values = colour)
plot_cell_trajectory(HSMM_primary, color_by = "SLC12A5_macro") +
  facet_wrap(~celltype, nrow = 1)+ scale_color_manual(values = colour)

plot_cell_trajectory(HSMM_primary,color_by = 'celltype')+scale_color_npg()
HSMM_primary_root_state_4 <- orderCells(HSMM_primary,root_state = 4) 

library(ggpubr)
df <- pData(HSMM_primary)

View(df)
ggplot(df,aes(Pseudotime,colour=SLC12A5_macro,fill=SLC12A5_macro))+
  geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic2()+
  theme(text = element_text(size = 20))

p <- plot_genes_jitter(HSMM_primary['SLC12A5',],color_by = 'State')
pData(HSMM_primary)$SLC12A5=log2(exprs(HSMM_primary)['SLC12A5',]+1)
plot_cell_trajectory(HSMM_primary,color_by = 'SLC12A5')+scale_color_gsea()
pData(HSMM_primary)$SPP1=log2(exprs(HSMM_primary)['SPP1',]+1)
plot_cell_trajectory(HSMM_primary,color_by = 'SPP1')+scale_color_gsea()



HSMM_all_SCT <- orderCells(HSMM_all_SCT,root_state = 9) 
pdf('20230508_HSMM_ALL_SCT_Pseudotime.pdf',height=5,width=5)
plot_cell_trajectory(HSMM_all_SCT,show_tree=T, color_by = "Pseudotime",show_branch_points=F)+scale_color_gsea()+ theme(text = element_text(size =20))
dev.off()
pdf('20230508_HSMM_ALL_SCT_celltype.pdf',height=5,width=5)
plot_cell_trajectory(HSMM_all_SCT,color_by = 'celltype',show_branch_points=F)+scale_color_npg()+ theme(text = element_text(size =20))
dev.off()

HSMM_primary <- readRDS('20230501_HSMM_all_SCT.rds')



macro <- readRDS("/data/shpcv2_kvcka/20230426_PDAC/Macro/20230501_macro_annotated.rds")
macro <- FindAllMarkers(macro,only.pos = T)
save(sce.markers,file = '20230508_macro.markers.Rdata')
library(patchwork)
sce.markers <- macro

marker_selected_1 <- sce.markers  %>% 
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::group_by(cluster) %>% 
  dplyr::slice_max(order_by = avg_log2FC, n = 10)



markers <- as.data.frame(marker_selected_1)
markerdata <- ScaleData(macro, features = as.character(unique(markers$gene)), assay = "RNA")

DoHeatmap(macro,
          features = as.character(unique(markers$gene)),
          group.by = "SLC12A5_macro",
          assay = 'RNA')
pdf('20230508_heatmap.pdf')
DoHeatmap(markerdata,label=T,
          features = as.character(unique(markers$gene)),
          group.by = "SLC12A5_macro",
          assay = 'RNA',
          group.colors = c("#00BFC4","#AB82FF","#00CD00","#C77CFF"))+
  scale_fill_gradientn(colors = c("white","grey","firebrick3"))+ theme(text = element_text(size =10))

dev.off()

color <- c('lightgrey','darkred')
pdf('20230510_VEGFA_FeaturePlot.pdf')
FeaturePlot(macro,order=T,reduction='tsne',features = c('VEGFA'),ncol = 1,cols=color,pt.size=1)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+ theme(text = element_text(size =30))
dev.off()

pdf('20230510_PIGF_FeaturePlot.pdf')
FeaturePlot(macro,order=T,reduction='tsne',features = c('PIGF'),ncol = 1,cols=color,pt.size=1)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+ theme(text = element_text(size =30))
dev.off()

pdf('20230510_CCL2_FeaturePlot.pdf')
FeaturePlot(macro,order=T,reduction='tsne',features = c('CCL2'),ncol = 1,cols=color,pt.size=1)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+ theme(text = element_text(size =30))
dev.off()

pdf('20230510_STAT3_FeaturePlot.pdf')
FeaturePlot(macro,order=T,reduction='tsne',features = c('STAT3'),ncol = 1,cols=color,pt.size=1)+
  theme(panel.border=element_rect(fill=NA,color='black',size=1,linetype='solid'))+ theme(text = element_text(size =30))
dev.off()
