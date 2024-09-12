setwd("/home/shpcv2_kvcka/data/20230426_PDAC")

load("/home/shpcv2_kvcka/data/20230426_PDAC/rawdata_list.Rdata")


obj.list <- rawdata_list
obj.list <- lapply(names(obj.list), function(x) {
  CreateSeuratObject(counts = obj.list[[x]],
                     project  = x,
                     min.cells = 3,     
                     min.features = 200)})  
lapply(obj.list,dim)
lapply(rawdata_list,dim)
names(obj.list) <- names(rawdata_list)
obj.list


obj.list$MET01@meta.data$group <- 'metastasis'
obj.list$MET02@meta.data$group <- 'metastasis'
obj.list$MET03@meta.data$group <- 'metastasis'
obj.list$MET04@meta.data$group <- 'metastasis'
obj.list$MET05@meta.data$group <- 'metastasis'
obj.list$MET06@meta.data$group <- 'metastasis'
obj.list$P01@meta.data$group <- 'primary'
obj.list$P02@meta.data$group <- 'primary'
obj.list$P03@meta.data$group <- 'primary'
obj.list$P04@meta.data$group <- 'primary'
obj.list$P05@meta.data$group <- 'primary'
obj.list$P06@meta.data$group <- 'primary'
obj.list$P07@meta.data$group <- 'primary'
obj.list$P08@meta.data$group <- 'primary'
obj.list$P09@meta.data$group <- 'primary'
obj.list$P10@meta.data$group <- 'primary'
obj.list$AdjNorm_TISSUE_1@meta.data$group <- 'adjnromal'
obj.list$AdjNorm_TISSUE_2@meta.data$group <- 'adjnromal'
obj.list$AdjNorm_TISSUE_3@meta.data$group <- 'adjnromal'
obj.list$PDAC_TISSUE_1@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_10@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_11A@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_11B@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_12@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_13@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_15@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_16@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_2@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_3@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_4@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_5@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_6@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_7@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_8@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_9@meta.data$group <- 'primary'
obj.list$PDAC_TISSUE_14@meta.data$group <- 'primary'
obj.list$Liver@meta.data$group <- 'metastasis'
obj.list$Lung@meta.data$group <- 'metastasis'
obj.list$Vaginal@meta.data$group <- 'metastasis'
obj.list$Peritoneal@meta.data$group <- 'metastasis'

obj.list$MET01@meta.data$GSM <- 'GSM4679542'
obj.list$MET02@meta.data$GSM <- 'GSM4679543'
obj.list$MET03@meta.data$GSM <- 'GSM4679544'
obj.list$MET04@meta.data$GSM <- 'GSM4679545'
obj.list$MET05@meta.data$GSM <- 'GSM4679546'
obj.list$MET06@meta.data$GSM <- 'GSM4679547'
obj.list$P01@meta.data$GSM <- 'GSM4679532'
obj.list$P02@meta.data$GSM <- 'GSM4679533'
obj.list$P03@meta.data$GSM <- 'GSM4679534'
obj.list$P04@meta.data$GSM <- 'GSM4679535'
obj.list$P05@meta.data$GSM <- 'GSM4679536'
obj.list$P06@meta.data$GSM <- 'GSM4679537'
obj.list$P07@meta.data$GSM <- 'GSM4679538'
obj.list$P08@meta.data$GSM <- 'GSM4679539'
obj.list$P09@meta.data$GSM <- 'GSM4679540'
obj.list$P10@meta.data$GSM <- 'GSM4679541'
obj.list$AdjNorm_TISSUE_1@meta.data$GSM <- 'GSM4710706'
obj.list$AdjNorm_TISSUE_2@meta.data$GSM <- 'GSM4710707'
obj.list$AdjNorm_TISSUE_3@meta.data$GSM <- 'GSM4710708'
obj.list$PDAC_TISSUE_1@meta.data$GSM <- 'GSM4710689'
obj.list$PDAC_TISSUE_10@meta.data$GSM <- 'GSM4710698'
obj.list$PDAC_TISSUE_11A@meta.data$GSM <- 'GSM4710699'
obj.list$PDAC_TISSUE_11B@meta.data$GSM <- 'GSM4710700'
obj.list$PDAC_TISSUE_12@meta.data$GSM <- 'GSM4710701'
obj.list$PDAC_TISSUE_13@meta.data$GSM <- 'GSM4710702'
obj.list$PDAC_TISSUE_15@meta.data$GSM <- 'GSM4710704'
obj.list$PDAC_TISSUE_16@meta.data$GSM <- 'GSM4710705'
obj.list$PDAC_TISSUE_2@meta.data$GSM <- 'GSM4710690'
obj.list$PDAC_TISSUE_3@meta.data$GSM <- 'GSM4710691'
obj.list$PDAC_TISSUE_4@meta.data$GSM <- 'GSM4710692'
obj.list$PDAC_TISSUE_5@meta.data$GSM <- 'GSM4710693'
obj.list$PDAC_TISSUE_6@meta.data$GSM <- 'GSM4710694'
obj.list$PDAC_TISSUE_7@meta.data$GSM <- 'GSM4710695'
obj.list$PDAC_TISSUE_8@meta.data$GSM <- 'GSM4710696'
obj.list$PDAC_TISSUE_9@meta.data$GSM <- 'GSM4710697'
obj.list$PDAC_TISSUE_14@meta.data$GSM <- 'GSM4710703'
obj.list$Liver@meta.data$GSM <- 'GSM4730266'
obj.list$Lung@meta.data$GSM <- 'GSM4730267'
obj.list$Vaginal@meta.data$GSM <- 'GSM4730265'
obj.list$Peritoneal@meta.data$GSM <- 'GSM4730268'


for (x in names(obj.list)){
  obj.list[[x]][["percent.MT"]] <- PercentageFeatureSet(obj.list[[x]], pattern = "^MT-")
  obj.list[[x]][["percent.RP"]] <- PercentageFeatureSet(obj.list[[x]], pattern = "^RP[SL]")
  obj.list[[x]][["percent.HB"]] <- PercentageFeatureSet(obj.list[[x]], pattern = "^Hb[^(p)]")
}

qc_feature <- c("nFeature_RNA", "nCount_RNA", "percent.HB", "percent.MT",  "percent.RP")

for(x in names(obj.list)){
  pdf(file=paste0("1_",x,"_quality_control.pdf"),width = 15,height=7)
  print(VlnPlot(obj.list[[x]], features = qc_feature, ncol = 5, pt.size = 0.5))
  dev.off()  
}

saveRDS(obj.list,'20230427_largelist.rds')

setwd("/home/shpcv2_kvcka/data/20230426_PDAC")
getwd()
obj.list <- readRDS("20230427_largelist.rds")
obj.list <- lapply(obj.list, function(x) {
  subset(x, subset = nFeature_RNA > 200 & 
           nFeature_RNA < 7500 &
           percent.MT< 20)})
obj.list





