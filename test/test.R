## ------------------------------------------------
## Usage example
## ------------------------------------------------
library(Seurat)
library(Matrix)
library(FNN)
library(igraph)
library(reticulate)
library(leiden)
library(dplyr)
library(ggplot2)
set.seed(42)      

source('FindSpaCC.R')
source('SpaCC.R')
source('SpaCCPlot.R')

setwd('E:\\stTNBC\\01')
seo <- CreateSeuratObject(counts = Read10X('./'), project = 'Posterior1', assay = 'Spatial')
seo$slice <- 1
seo$region <- 'Posterior1' 

imgpath <- "./Spatial/"
img <- Seurat::Read10X_Image(image.dir = imgpath)
Seurat::DefaultAssay(object = img) <- 'Spatial'
img <- img[colnames(x = seo)]
seo[['image']] <- img
seo
seo <- Seurat::SCTransform(seo, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)
seo <- Seurat::RunPCA(seo, assay = "SCT", verbose = FALSE)
seo <- Seurat::FindNeighbors(seo, reduction = "pca", dims = 1:30)
seo <- Seurat::FindClusters(seo, verbose = FALSE)
seo <- Seurat::RunUMAP(seo, reduction = "pca", dims = 1:30)
save(seo,file = '01.robj')

rm(list = ls())
gc()

load('01.robj')
seo
colnames(seo@meta.data)

##Spatial coor
spatial_coords <-seo@images$image@coordinates[,c('imagerow','imagecol')]
head(spatial_coords)

seo <- SpaCC(
  seo = seo,
  spatial_coords = spatial_coords
)

pdf('01_SpaCC.pdf')
SpaCCPlot(seo,feature = 'SpaCC')
dev.off()

library(fossil)
ari <- fossil::adj.rand.index(seo$SpaCC, seo$seurat_clusters)
ari

library(profvis)
# 使用 profvis 分析性能
profvis({
  # 这里是你要测试的代码
  seo <- SpaCC(
    seo = seo,
    spatial_coords = spatial_coords
  )
})
