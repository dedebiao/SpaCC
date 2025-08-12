## ------------------------------------------------
## Usage example
## ------------------------------------------------
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

library(profvis)
# 使用 profvis 分析代码性能
profvis({
  # 这里是你要测试的代码
  seo <- SpaCC(
    seo = seo,
    spatial_coords = spatial_coords
  )
})
