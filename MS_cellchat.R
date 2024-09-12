rm(list=ls())
setwd("/home/shpcv2_kvcka/data/20230426_PDAC/MS_cellchat")
getwd()
library(CellChat)
library(patchwork)

options(stringsAsFactors = FALSE)
Idents(MS) <- "group"
MS <- subset(MS,idents = "primary")
Idents(MS) <- "SLC12A5_MS"
DimPlot(MS,reduction = 'tsne')

cellchat <- createCellChat(MS@assays$SCT@data)
MS$cellType <- Idents(MS) 
meta <- data.frame(cellType = MS$cellType, row.names =  Cells(MS))
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
pdf('20230507_MS_cellchat_all.pdf',width =6,height =6)
netVisual_circle(cellchat_all@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_all@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


save(cellchat_all,file = '20230505_cellchat_all_result.Rdata')


rm(list=ls())
gc()
setwd("/home/shpcv2_kvcka/data/20230426_PDAC/MS_cellchat")
getwd()
load("/data/shpcv2_kvcka/20230426_PDAC/MS_cellchat/20230507_MS_cellchat_all_result.Rdata")
pdf('20230507_MS_bubbleplot_all.pdf',height=39,width=20)
netVisual_bubble(cellchat_all, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)
dev.off()
pdf('20230507_MS_bubbleplot_all_source_SLC12A5.pdf',height=18,width=10)
netVisual_bubble(cellchat_all, sources.use = c('mM2-TAMs',"bM2-TAMs"), targets.use = NULL, remove.isolate = FALSE)+theme(text = element_text(size = 18))
dev.off()
pdf('20230507_MS_bubbleplot_all_target_SLC12A5.pdf',height=24,width=10)
netVisual_bubble(cellchat_all, sources.use = NULL, targets.use = c('mM2-TAMs',"bM2-TAMs"), remove.isolate = FALSE)+theme(text = element_text(size = 18))
dev.off()
cellchat_all@netP$pathways
netVisual_aggregate(cellchat_all, signaling = 'PTPRM', layout = "circle")
plotGeneExpression(cellchat_all, signaling = "PDGF")
cellchat_all <- netAnalysis_computeCentrality(cellchat_all, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat_all, signaling = 'SPP1', width = 8, height = 2.5, font.size = 10)

pdf('20230507_netAnalysis_signalingRole_scatter_cellchat_all.pdf')
netAnalysis_signalingRole_scatter(cellchat_all,label.size = 7)+theme(text = element_text(size = 18))
dev.off()

cellchat_all <- computeNetSimilarity(cellchat_all, type = "functional")
cellchat_all <- netEmbedding(cellchat_all, type = "functional",umap.method = 'uwot')

future::plan("multisession", workers = 4) # do parallel
cellchat_all <- netClustering(cellchat_all, type = "functional",do.parallel=F)

pdf('20230507_netVisual_embedding_cellchat_all.pdf')
netVisual_embedding(cellchat_all, type = "functional", label.size =4.2)+theme(text = element_text(size = 18))
dev.off()
save(cellchat_all,file='20230507_MS_cellchat_all_final.Rdata')


cellchat_all@netP$pathways
pathway.list <- c('VCAM','TNF','SPP1',
                  'OSM','IL16','GAS')
#'THBS','RETN
'F1' %in% cellchat_all@netP$pathways
plotGeneExpression(cellchat_all, signaling = "COLLAGEN")

for (i in pathway.list){
  dir.create(paste0('cellchatplot/',i))
  pdf(paste0('cellchatplot/',i,'/20230507_',i,'_networkplot.pdf'),height=5,width = 5)
  netAnalysis_signalingRole_network(cellchat_all, signaling = i, width = 8, height = 2.5, font.size = 10,font.size.title = 13)
  dev.off()
  
  pdf(paste0('cellchatplot/',i,'/20230507_',i,'_Heatmapplot.pdf'),height = 4,width=4)
  p1 <- netVisual_heatmap(cellchat_all, signaling = i, color.heatmap = "Reds")
  print(p1)
  dev.off()
  LR.show <- extractEnrichedLR(cellchat_all, signaling = i, geneLR.return = FALSE) # show one ligand-receptor pair
  pdf(paste0('cellchatplot/',i,'/20230507_',i,'_circleplot.pdf'),height=30,width=30)
  netVisual_individual(cellchat_all, signaling = i, pairLR.use = LR.show, layout = "circle",nCol=4)
  dev.off()
  pdf(paste0('cellchatplot/',i,'/20230507_',i,'_chordplot.pdf'))
  p2 <- netVisual_individual(cellchat_all, signaling = i, pairLR.use = LR.show, layout = "chord",nCol=2)
  print(p2)
  dev.off()
  pdf(paste0('cellchatplot/',i,'/20230507_',i,'_contributionplot.pdf'))
  p3 <- netAnalysis_contribution(cellchat_all, signaling = i)
  print(p3)
  dev.off()
  pdf(paste0('cellchatplot/',i,'/20230507_',i,'_Vlnplot.pdf'),height = 4,width=4)
  p4 <- plotGeneExpression(cellchat_all, signaling = i)
  print(p4)
  dev.off()
  
}
