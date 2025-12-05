library(Seurat)
library(harmony)
library(reticulate)
count<-read.table("data/singlecell/21762584/scmatrix.tsv.gz",row.names = 1,header = T)
meta<-read.table("data/singlecell/21762584/scmetadata10.txt.gz",row.names = 1,sep = "\t")
meta<-meta[colnames(count),]
sce.all<-CreateSeuratObject(count,)

VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,raster=FALSE)
plot1 <- FeatureScatter(sce.all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

sce.all <- subset(sce.all, subset = nFeature_RNA > 400 & nFeature_RNA < 4000 & 
                    percent.mt < 10&nCount_RNA<25000&nCount_RNA>5000)

sce.all <- NormalizeData(sce.all, normalization.method = "LogNormalize", scale.factor = 10000)
sce.all <- FindVariableFeatures(sce.all, selection.method = "vst", nfeatures = 5000)

sce.all <- ScaleData(sce.all)
sce.all <- RunPCA(sce.all, features = VariableFeatures(object = sce.all))

scRNA_harmony <- RunHarmony(sce.all,reduction = "pca",group.by.vars = "group",reduction.save = "harmony")
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:15,reduction.name = "umap")
DimPlot(scRNA_harmony, reduction = "umap",group.by = "group",raster=FALSE)
              
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:15) %>%  
  FindClusters(resolution = 0.3) 
DimPlot(scRNA_harmony, reduction = "umap",group.by = "group",raster=FALSE)
DimPlot(scRNA_harmony, reduction = "umap",group.by = "seurat_clusters",raster=FALSE,label = T)
scRNA_harmony<-JoinLayers(scRNA_harmony)
write.csv(t(scRNA_harmony@assays$RNA$counts),"data/singlecell/matrix.csv")

celltypist<-import("celltypist")
models<-celltypist$models
models$models_path
models$download_models(model = 'Adult_COVID19_PBMC.pkl')
model = models$Model$load(model = 'Adult_COVID19_PBMC.pkl')
predictions = celltypist$annotate("data/singlecell/matrix.csv", model = model)
scRNA_harmony$celltype<-predictions$predicted_labels

DimPlot(scRNA_harmony, reduction = "umap",group.by = "celltype",raster=FALSE,label = T)
