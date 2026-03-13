library(Seurat)
library(harmony)
library(ggplot2)
library(scop)
library(DoubletFinder)
library(clustree)
meta.data<-read.csv("data/SraRunTable.csv",header = T)
sample<-list.files("data/singlecell/",pattern = "SRR")
seu_list<-list()
for(i in 1:length(sample)){
  path<-paste0("data/singlecell/",sample[i])
  count<-Read10X(path)
  seu<-CreateSeuratObject(count,project = sample[i],min.cells = 20,min.features = 200)
  seu$age<-meta.data$AGE[i]
  seu$sex<-meta.data$sex[i]
  seu$group<-ifelse(grepl("Control",meta.data$Library.Name[i]),"Control","Hypertension")
  seu$orig.ident<-meta.data$Library.Name[i]
  seu_list[[i]]<-seu
}

seu_list<-lapply(seu_list,function(x){
  x[["percent.mt"]] = PercentageFeatureSet(x, pattern = "^MT-")
  pbmc<-subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                 percent.mt < 10&nCount_RNA<15000&nCount_RNA>1000)
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc<- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  scale.genes <-  VariableFeatures(pbmc)
  pbmc <- ScaleData(pbmc, features = scale.genes)
  pbmc<- RunPCA(pbmc, features = VariableFeatures(pbmc))
  pbmc<-RunUMAP(pbmc, dims = 1:20,reduction = 'pca')
  pbmc <- FindNeighbors(pbmc, dims = 1:20, reduction = "pca")
  pbmc <- FindClusters(pbmc,resolution = 0.6)
  sweep.res.list <- paramSweep(JoinLayers(pbmc), PCs = 1:20, sct = FALSE,num.cores = 8)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  homotypic<-modelHomotypic(pbmc@meta.data$seurat_clusters)
  Doubletrate <- 0.05
  nExp_poi <- round(Doubletrate*ncol(pbmc))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic))
  pbmc_douletfinder <- doubletFinder(pbmc, PCs = 1:20, pN = 0.25, pK = mpK, 
                                     nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)
  pbmc_douletfinder
})

seu_list<-lapply(seu_list,function(x){
  x[,x@meta.data[,11]=="Singlet"]
})
saveRDS(seu_list,"data/seu_list.rds")
seu_list<-readRDS("data/singlecell/seu_list.rds")



pbmc<-merge(seu_list[[1]],seu_list[-1])

Idents(pbmc)<-"orig.ident"

#数据过滤----
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

seu_list[[1]]<-subset(seu_list[[1]], subset = nFeature_RNA < 3300 &nCount_RNA<11000)
seu_list[[2]]<-subset(seu_list[[2]], subset = nFeature_RNA < 3100 &nCount_RNA<10000)
seu_list[[3]]<-subset(seu_list[[3]], subset = nFeature_RNA < 3000 &nCount_RNA<9500)
seu_list[[4]]<-subset(seu_list[[4]], subset = nFeature_RNA < 3000 &nCount_RNA<9000)
seu_list[[5]]<-subset(seu_list[[5]], subset = nFeature_RNA < 2950 &nCount_RNA<9000)
seu_list[[6]]<-subset(seu_list[[6]], subset = nFeature_RNA < 3500 &nCount_RNA<10000)
#seu_list[[7]]<-subset(seu_list[[7]], subset = nFeature_RNA < 3500 &nCount_RNA<12000)
seu_list[[8]]<-subset(seu_list[[8]], subset = nFeature_RNA < 3500 &nCount_RNA<12000)
seu_list[[9]]<-subset(seu_list[[9]], subset = nFeature_RNA < 3300 &nCount_RNA<11000)
#seu_list[[10]]<-subset(seu_list[[10]], subset = nFeature_RNA < 3500 &nCount_RNA<12000)
seu_list[[11]]<-subset(seu_list[[11]], subset = nFeature_RNA < 3000 &nCount_RNA<10000)
#seu_list[[12]]<-subset(seu_list[[12]], subset = nFeature_RNA < 3500 &nCount_RNA<12000)
seu_list[[13]]<-subset(seu_list[[13]], subset = nFeature_RNA < 3300 &nCount_RNA<12000)
#seu_list[[14]]<-subset(seu_list[[14]], subset = nFeature_RNA < 3500 &nCount_RNA<12000)
seu_list[[15]]<-subset(seu_list[[15]], subset = nFeature_RNA < 3500 &nCount_RNA<11000)
seu_list[[16]]<-subset(seu_list[[16]], subset = nFeature_RNA < 3500 &nCount_RNA<12000)
seu_list[[17]]<-subset(seu_list[[17]], subset = nFeature_RNA < 3500 &nCount_RNA<11500)
seu_list[[18]]<-subset(seu_list[[18]], subset = nFeature_RNA < 3500 &nCount_RNA<12000)
seu_list[[19]]<-subset(seu_list[[19]], subset = nFeature_RNA < 3100 &nCount_RNA<11500)
seu_list[[20]]<-subset(seu_list[[20]], subset = nFeature_RNA < 3800 &nCount_RNA<12000)
seu_list[[21]]<-subset(seu_list[[21]], subset = nFeature_RNA < 2500 &nCount_RNA<8000)
seu_list[[22]]<-subset(seu_list[[22]], subset = nFeature_RNA < 3000 &nCount_RNA<9000)
seu_list[[23]]<-subset(seu_list[[23]], subset = nFeature_RNA < 2950 &nCount_RNA<9000)
seu_list[[24]]<-subset(seu_list[[24]], subset = nFeature_RNA < 2900 &nCount_RNA<9000)
seu_list[[25]]<-subset(seu_list[[25]], subset = nFeature_RNA < 2700 &nCount_RNA<7000)

pbmc<-merge(seu_list[[1]],seu_list[-1])
pbmc<-subset(pbmc,subset = nFeature_RNA>500)
Idents(pbmc)<-"orig.ident"
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc<- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
scale.genes <-  VariableFeatures(pbmc)
pbmc <- ScaleData(pbmc, features = scale.genes)
pbmc<- RunPCA(pbmc, features = VariableFeatures(pbmc))
DimPlot(pbmc, reduction = "pca", group.by = "orig.ident")

#去批次----
scRNA_harmony <- RunHarmony(pbmc, group.by.vars = "orig.ident")
DimPlot(scRNA_harmony, reduction = "harmony", group.by = "orig.ident",raster = F)

ElbowPlot(scRNA_harmony,reduction = 'harmony',ndims = 50)
scRNA_harmony <- RunUMAP(scRNA_harmony, dims = 1:30,reduction = 'harmony')
#scRNA_harmony<-RunTSNE(scRNA_harmony,dims = 1:30,reduction = "harmony")
pbmc <- FindNeighbors(scRNA_harmony, dims = 1:30, reduction = "harmony")
pbmc <- FindClusters(pbmc,resolution = seq(0.4,0.8,0.1))
clustree(pbmc@meta.data, prefix ="RNA_snn_res.")
DimPlot(pbmc,reduction = "umap",group.by = "RNA_snn_res.0.6",raster = F,label = T)

pbmc<-JoinLayers(pbmc)
Matrix::writeMM(t(pbmc@assays$RNA$counts),"data/singlecell/matrix.mtx")
write.table(row.names(pbmc),"data/singlecell/genefile.tsv",sep = "\t",col.names = F,row.names = F)
write.table(colnames(pbmc),"data/singlecell/cellfile.tsv",sep = "\t",col.names = F,row.names = F)

celltypist<-reticulate::import("celltypist")
models<-celltypist$models
models$models_path
models$download_models(model = 'Adult_COVID19_PBMC.pkl')
model = models$Model$load(model = 'Adult_COVID19_PBMC.pkl')
predictions = celltypist$annotate("data/singlecell/matrix.mtx", model = model,
                                  gene_file="data/singlecell/genefile.tsv",
                                  cell_file="data/singlecell/cellfile.tsv")
pbmc$celltype<-predictions$predicted_labels

DimPlot(pbmc, reduction = "umap",group.by = "celltype",raster=FALSE)
ggsave("result/celltype.pdf",width=7,height = 5)

saveRDS(pbmc,"data/pbmc.rds")

pbmc<-readRDS("data/pbmc.rds")
trace(scop:::GetAssayData5.Assay5,edit = T)
pbmc <- standard_scop(srt = pbmc)

pancreas_sub <- RunPAGA(
  srt = pbmc,
  group_by = "celltype",
  linear_reduction = "pca",
  nonlinear_reduction = "umap"
)
PAGAPlot(
  srt = pbmc,
  reduction = "umap",
  label = TRUE,
  label_insitu = TRUE,
  label_repel = TRUE
)
library(Seurat)
DefaultAssay(pbmc) <- "RNA"

library(scCustomize)
seurat_obj_v4 <- Convert_Assay(pbmc, convert_to = "V3",assay = "RNA")
SaveH5Seurat(seurat_obj_v4, filename = "data.h5Seurat",overwrite = T)

Convert("data.h5Seurat", dest = "h5ad", overwrite = TRUE)

cellann<-data.frame(cell_id=colnames(pbmc),
                    cell_annotation=pbmc$celltype)
write.table(cellann,"data/singlecell/cell_annotation.tsv",sep = "\t",quote = F,row.names = F)
