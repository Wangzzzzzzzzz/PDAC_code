setwd("/home/shpcv2_kvcka/data/20230426_PDAC/Immune")
getwd()
obj.combined <-readRDS('20230430_Immune_annotated.rds')

#DefaultAssay(obj.combined) <- "SCT"
M2.diff <- FindMarkers(obj.combined, ident.1 = "SLC12A5_high_M2", ident.2 = "SLC12A5_low_M2", logfc.threshold = 0,verbose = T)
write.csv(file="20230430_M2.diff.csv",M2.diff)

# Idents(obj.combined) <- 'celltype'
# Mac <- subset(obj.combined, idents = "Mac")
# Idents(Mac) <- "group"
# avg.Mac <- log1p(AverageExpression(Mac, verbose = FALSE)$RNA)#log1p(x) computes log(1+x) accurately also for |x| << 1.
# avg.Mac <- data.frame(avg.Mac ,gene=rownames(avg.Mac))

low<-floor(range(M2.diff$avg_log2FC)[1]) 
high<-ceiling(range(M2.diff$avg_log2FC)[2]) 
pdf('20230430_DEG_M2.pdf',width = 10,height=10)
print(EnhancedVolcano(M2.diff,
                      title = 'SLC12A5_high_M2 versus SLC12A5_low_M2',
                      lab = rownames(M2.diff),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      FCcutoff = 0.5,
                      xlim = c(low, high))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()))
dev.off()

## KEGG pathway analysis
R.utils::setOption("clusterProfiler.download.method",'auto')

rt <- M2.diff
rt$id <- rownames(rt)
rt <- subset(rt,subset=p_val<0.05)
rt <- subset(rt,subset=abs(avg_log2FC)>0.25)
rt <- rt[,c(-1,-3,-4,-5)]
genes=as.vector(rt[,2])  
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)   
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)

gene_up <- subset(out,subset=avg_log2FC>0)
gene_down <- subset(out,subset=avg_log2FC<0)
R.utils::setOption("clusterProfiler.download.method",'auto')
kk.up <- enrichKEGG(gene         = gene_up$entrezID,
                    organism     = 'hsa',
                    #universe     = gene_all,
                    pvalueCutoff = 0.05,
                    qvalueCutoff =0.05)
head(kk.up)[,1:6]

kk.down <- enrichKEGG(gene         =  gene_down$entrezID,
                      organism     = 'hsa',
                      #universe     = gene_all,
                      pvalueCutoff = 0.05,
                      qvalueCutoff =0.05)
head(kk.down)[,1:6]

kk.up = setReadable(kk.up,
                    OrgDb = "org.Hs.eg.db",
                    keyType = "ENTREZID")

kk.down = setReadable(kk.down,
                      OrgDb = "org.Hs.eg.db",
                      keyType = "ENTREZID")

kegg_down_dt <- as.data.frame(kk.down)
kegg_up_dt <- as.data.frame(kk.up)
down_kegg<-kegg_down_dt[kegg_down_dt$pvalue<0.05,];down_kegg$group=-1
up_kegg<-kegg_up_dt[kegg_up_dt$pvalue<0.05,];up_kegg$group=1

write.csv(file="20230430_SLC12A5HighM2.vs.SLC12A5LowM2.diff_KEGG_up.csv",kegg_up_dt)
write.csv(file="20230430_SLC12A5HighM2.vs.SLC12A5LowM2.diff_KEGG_down.csv",kegg_down_dt)



library(ggplot2)
library(ggsci)
dat=rbind(up_kegg,down_kegg)
for (i in 1:nrow(dat)){
  print(i)
  dat$Description[i] <- substr(as.character(dat$Description[i]),1,(nchar(as.character(dat$Description[i]))-23))
}
dat$pvalue = -log10(dat$pvalue)
dat[dat$group==-1,]$pvalue <- dat[dat$group==-1,]$pvalue*-1
dat$description <- dat$Description
dat[duplicated(dat$Description),"Description"] <- paste0(dat[duplicated(dat$Description),"Description"],".1") # 调整y轴
dat$Description <- factor(dat$Description,levels = dat[order(dat$pvalue,decreasing = T),"Description"]) # 设置因子
name.y <- substr(dat$Description,1,1)
names(name.y) <- dat$description

ggplot(dat,aes(x=reorder(Description,order(pvalue, decreasing = T)),y=pvalue, fill=factor(group)))+geom_bar(stat="identity")+coord_flip( )+
  scale_x_discrete(
    "Description",
    labels = names(name.y))+
  scale_fill_d3(breaks=c('-1','1'),labels=c('down','up'),name='gene group')+
  theme(axis.text=element_text(face = "bold",size = 15),
        axis.title = element_text(face = 'bold',size = 15),
        plot.title = element_text(size = 20,hjust = 0.3),
        legend.position = "top",
        panel.grid = element_line(colour = 'white'))

pdf(file="20230430_SLC12A5HighM2.vs.SLC12A5LowM2.diff_KEGG.pdf",width=10,height = 15)
last_plot()+theme(legend.key.size=unit(0.2,'inches'),legend.text = element_text(size = 15),legend.title=element_text(size=16))#调整图例大小
dev.off()

rt <- M2.diff
rt$id <- rownames(rt)
rt <- subset(rt,subset=p_val<0.1)
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
gseaplot(kk_gse, geneSetID = rownames(kk_gse[1,]))
write.csv(file="20230430_SLC12A5HighM2.vs.SLC12A5LowM2.diff_GSEA_KEGG_result.csv",data.frame(kk_gse))
pdf(file="20230430_SLC12A5HighM2.vs.SLC12A5LowM2.diff_GSEA_KEGG.pdf",height = 7)
clusterProfiler::dotplot(kk_gse,showCategory=20,label_format=80) 
dev.off()

devtools::install_github("junjunlab/GseaVis")
library(GseaVis)
pdf(file="macro/20230106_GSEA_KEGG/202301010_GSEA_hsa04066.pdf",width=12)
p1 <- gseaNb(object = kk_gse,
             geneSetID = 'hsa04066',
             newGsea = F)

dev.off()

