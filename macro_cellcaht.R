rm(list=ls())
setwd("/home/shpcv2_kvcka/data/20230426_PDAC/Immune")
getwd()
# Immune <- readRDS('20230503_Immune_annotated_with_SLC12A5_macro.rds')


table(Immune@active.ident)
paste0(sprintf("'%s'", levels(Immune)), collapse = ",")
new.cluster.ids <- c( 'mM2-TAMs','M1-TAMs','Cycle',
                      'bM2-TAMs','CD4_T','DC','Treg',
                      'unknown','CD8_T','B_cells','NK',
                      'Neutrophils','Plasma','Mast' )

names(new.cluster.ids) <- levels(Immune)
Immune <- RenameIdents(Immune, new.cluster.ids)
table(Immune@active.ident)
Immune$mbcelltype <- Immune@active.ident
saveRDS(Immune,'20230505_Immune_mbcelltype_annotated.rds')



rm(list=ls())
setwd("/home/shpcv2_kvcka/data/20230426_PDAC/Immune")
getwd()
Immune <- readRDS('20230505_Immune_mbcelltype_annotated.rds')

setwd("/home/shpcv2_kvcka/data/20230426_PDAC/Immune_cellchat")
getwd()

library(CellChat)
library(patchwork)

options(stringsAsFactors = FALSE)
Idents(Immune) <- "group"
Immune <- subset(Immune,idents = "primary")
Idents(Immune) <- "mbcelltype"
Immune <- subset(Immune,subset=mbcelltype!='unknown')
Immune <- subset(Immune,subset=mbcelltype!='Cycle')
DimPlot(Immune,reduction = 'tsne')

cellchat <- createCellChat(Immune@assays$SCT@data)
Immune$cellType <- Idents(Immune) 
meta <- data.frame(cellType = Immune$cellType, row.names =  Cells(Immune))
cellchat <- addMeta(cellchat, meta = meta, meta.name = "cellType")

cellchat <- setIdent(cellchat, ident.use = "cellType") 
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

cellchat_all <- cellchat
pdf('20230505_cellchat_all.pdf',width =6,height =6)
netVisual_circle(cellchat_all@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_all@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


save(cellchat_all,file = '20230505_cellchat_all_result.Rdata')

cellchat <- createCellChat(Immune@assays$SCT@data)
Immune$cellType <- Idents(Immune) 
meta <- data.frame(cellType = Immune$cellType, row.names =  Cells(Immune))
cellchat <- addMeta(cellchat, meta = meta, meta.name = "cellType")

cellchat <- setIdent(cellchat, ident.use = "cellType") 
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
cellchat_ss <- cellchat
pdf('20230504_cellchat_all.pdf',width =6,height =6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


cellchat <- createCellChat(Immune@assays$SCT@data)
Immune$cellType <- Idents(Immune)
meta <- data.frame(cellType = Immune$cellType, row.names =  Cells(Immune))
cellchat <- addMeta(cellchat, meta = meta, meta.name = "cellType")

cellchat <- setIdent(cellchat, ident.use = "cellType") 
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

cellchat_cc <- cellchat
pdf('20230504_cellchat_cc.pdf',width =6,height =6)
netVisual_circle(cellchat_cc@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_cc@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()



cellchat <- createCellChat(Immune@assays$SCT@data)
Immune$cellType <- Idents(Immune) 
meta <- data.frame(cellType = Immune$cellType, row.names =  Cells(Immune))
cellchat <- addMeta(cellchat, meta = meta, meta.name = "cellType")

cellchat <- setIdent(cellchat, ident.use = "cellType") 
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)

cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

cellchat_er <- cellchat
pdf('20230504_cellchat_er.pdf',width =6,height =6)
netVisual_circle(cellchat_er@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_er@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()



save(cellchat_ss,cellchat_cc,cellchat_er,file = '20230504_cellchat_result.Rdata')


rm(list=ls())
gc()
setwd("/home/shpcv2_kvcka/data/20230426_PDAC/Immune_cellchat")
getwd()
load('20230505_cellchat_all_result.Rdata')

pdf('20230505_bubbleplot_all.pdf',height=15,width=20)
netVisual_bubble(cellchat_all, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)
dev.off()
pdf('20230505_bubbleplot_all_source_SLC12A5.pdf',height=18,width=15)
netVisual_bubble(cellchat_all, sources.use = c('mM2-TAMs',"bM2-TAMs"), targets.use = NULL, remove.isolate = FALSE)+theme(text = element_text(size = 18))
dev.off()
cellchat_all@netP$pathways
netVisual_aggregate(cellchat_all, signaling = 'SPP1', layout = "circle")
plotGeneExpression(cellchat_all, signaling = "SPP1")
cellchat_all <- netAnalysis_computeCentrality(cellchat_all, slot.name = "netP")
pdf('')
netAnalysis_signalingRole_network(cellchat_all, signaling = 'SPP1', width = 8, height = 2.5, font.size = 10)

pdf('20230505_netAnalysis_signalingRole_scatter_cellchat_all.pdf')
netAnalysis_signalingRole_scatter(cellchat_all,label.size = 7)+theme(text = element_text(size = 18))
dev.off()

cellchat_all <- computeNetSimilarity(cellchat_all, type = "functional")
cellchat_all <- netEmbedding(cellchat_all, type = "functional",umap.method = 'uwot')

future::plan("multisession", workers = 4) 
cellchat_all <- netClustering(cellchat_all, type = "functional",do.parallel=F)

pdf('20230505_netVisual_embedding_cellchat_all.pdf')
netVisual_embedding(cellchat_all, type = "functional", label.size = 3.5)
dev.off()
save(cellchat_all,file='20230505_cellchat_all_final.Rdata')



pairLR.CADM <- extractEnrichedLR(cellchat_all, signaling = 'CADM', geneLR.return = FALSE)
LR.show <- pairLR.CADM[1,] 

vertex.receiver = seq(1,4) 
netVisual_individual(cellchat_all, signaling = 'CADM',  pairLR.use = LR.show, vertex.receiver = vertex.receiver,layout = "hierarchy")

pathways.show <- c('CADM') 

vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat_all, signaling = pathways.show,  vertex.receiver = vertex.receiver,layout = "hierarchy")#layout = c("circle", "hierarchy", "chord", "spatial")


cellchat_all@netP$pathways
pathway.list <- c('SPP1','MIF','GAS','FN1',
                  'CCL','CADM','COMPLEMENT','ANNEXIN','ALCAM')


for (i in pathway.list){
  dir.create(paste0('cellchatplot/',i))
  pdf(paste0('cellchatplot/',i,'/20230506_',i,'_networkplot.pdf'),height=5,width = 5)
  netAnalysis_signalingRole_network(cellchat_all, signaling = i, width = 8, height = 2.5, font.size = 10,font.size.title = 13)
  dev.off()
  
  pdf(paste0('cellchatplot/',i,'/20230506_',i,'_Heatmapplot.pdf'),height = 4,width=4)
  p1 <- netVisual_heatmap(cellchat_all, signaling = i, color.heatmap = "Reds")
  print(p1)
  dev.off()
  LR.show <- extractEnrichedLR(cellchat_all, signaling = i, geneLR.return = FALSE) 
  pdf(paste0('cellchatplot/',i,'/20230506_',i,'_circleplot.pdf'))
  netVisual_individual(cellchat_all, signaling = i, pairLR.use = LR.show, layout = "circle",nCol=2)
  dev.off()
  pdf(paste0('cellchatplot/',i,'/20230506_',i,'_chordplot.pdf'))
  p2 <- netVisual_individual(cellchat_all, signaling = i, pairLR.use = LR.show, layout = "chord",nCol=2)
  print(p2)
  dev.off()
  pdf(paste0('cellchatplot/',i,'/20230506_',i,'_contributionplot.pdf'))
  p3 <- netAnalysis_contribution(cellchat_all, signaling = i)
  print(p3)
  dev.off()
  pdf(paste0('cellchatplot/',i,'/20230506_',i,'_Vlnplot.pdf'),height = 4,width=4)
  p4 <- plotGeneExpression(cellchat_all, signaling = i)
  print(p4)
  dev.off()
  
  }


netVisual_heatmap(cellchat_all, signaling = pathways.show, color.heatmap = "Reds")

pairLR.CXCL <- extractEnrichedLR(cellchat_all, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] 
pdf(paste0('cellchatplot/',pathways.show,'/20230506_',pathways.show,'_circleplot.pdf'),height=20,width = 20)
netVisual_individual(cellchat_all, signaling = pathways.show, pairLR.use = LR.show, layout = "circle",nCol=3)
dev.off()
