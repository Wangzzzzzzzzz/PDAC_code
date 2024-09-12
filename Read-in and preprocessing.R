# install.packages('devtools') 
# install.packages('remotes')
# install.packages('spatstat.explore') 
# remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
# install.packages('Seurat')
# if (!require("stringr", quietly = TRUE))
#   install.packages("stringr")
# if (!require("dplyr", quietly = TRUE))
#   install.packages("dplyr")
# sessionInfo('Seurat')
#install.packages('hdf5r')
# if (!require("Hmisc", quietly = TRUE))
#     install.packages("Hmisc")

library(devtools)
library(remotes)
library(stringr)
library(Seurat)
library(dplyr)
library(hdf5r)
library(Hmisc)
library(cowplot)
library(ggplot2)

setwd("/home/shpcv2_kvcka/data/20230426_PDAC/RawData/rawdata")
#1
path<-"GSE154778_RAW"
files <- list.files(path =path)    
fullpath<-paste(path,files,sep="/")
names(fullpath) <- lapply(files,function(x) {str_extract(x, '[\\w\\d]+')}) 
rawdata <- lapply(fullpath, function(x) {            
})
lapply(rawdata, dim)
#2
path<-"GSE155698_RAW"
files <- list.files(path =path)
files <- files[-11]
files <- files[-10]
fullpath<-paste(path,files,sep="/")
names(fullpath) <- lapply(files,function(x) {str_extract(x, '[\\w\\d]+')})   
rawdata2 <- lapply(fullpath, function(x) {             
  cells.table <- Read10X(x)
})
lapply(rawdata2, dim)
#3
path<-"GSE155698_RAW"
files <- list.files(path =path)
files <- files[10]
fullpath<-paste(path,files,'filtered_feature_bc_matrix.h5',sep="/")
names(fullpath) <- lapply(files,function(x) {str_extract(x, '[\\w\\d]+')})   
rawdata3 <- lapply(fullpath, function(x) {            
  cells.table <- Read10X_h5(x)
})
lapply(rawdata3, dim)
#4
path<-"GSE156405_Raw"
files <- list.files(path =path)
files <- files[-3]
fullpath<-paste(path,files,sep="/")
names(fullpath) <- lapply(files,function(x) {str_extract(x, '[\\w\\d]+')})  
rawdata4 <- lapply(fullpath, function(x) {            
  cells.table <- Read10X(x)
})
lapply(rawdata4, dim)
#5
path<-"GSE156405_Raw"
files <- list.files(path =path)
files <- files[3]
fullpath<-paste(path,files,sep="/")
names(fullpath) <- lapply(files,function(x) {str_extract(x, '[\\w\\d]+')})   
rawdata5 <- read.csv(fullpath,sep = ",",check.names = F, header = TRUE, row.names = 1)
rawdata5 <- t(rawdata5)
dim(rawdata5)
rawdata5 <- list(rawdata5)
names(rawdata5) <- c('Peritoneal')
lapply(rawdata5,dim)
setwd("/home/shpcv2_kvcka/data/20230426_PDAC")

rawdata_list <- c(rawdata,rawdata2,rawdata3,rawdata4,rawdata5)
lapply(rawdata_list,dim)

save(rawdata_list,file = 'rawdata_list.Rdata')
