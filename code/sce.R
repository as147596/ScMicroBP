library(Seurat)
library(harmony)
library(ggplot2)
library(scop)
meta.data<-read.csv("data/SraRunTable.csv",header = T)
sample<-list.files("data/singlecell/")
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

pbmc<-merge(seu_list[[1]],seu_list[-1])

pbmc[["percent.mt"]] = PercentageFeatureSet(pbmc, pattern = "^MT-")
Idents(pbmc)<-"orig.ident"

#数据过滤----
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & 
                    percent.mt < 10&nCount_RNA<15000&nCount_RNA>1000)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc<- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(pbmc)
pbmc <- ScaleData(pbmc, features = scale.genes)
pbmc<- RunPCA(pbmc, features = VariableFeatures(pbmc))
DimPlot(pbmc, reduction = "pca", group.by = "orig.ident")

#去批次----
scRNA_harmony <- RunHarmony(pbmc, group.by.vars = "orig.ident")
DimPlot(scRNA_harmony, reduction = "harmony", group.by = "orig.ident",raster = F)

ElbowPlot(scRNA_harmony,reduction = 'harmony',ndims = 50)
scRNA_harmony <- RunUMAP(scRNA_harmony, dims = 1:20,reduction = 'harmony')
scRNA_harmony<-RunTSNE(scRNA_harmony,dims = 1:20,reduction = "harmony")
pbmc <- FindNeighbors(scRNA_harmony, dims = 1:20, reduction = "harmony")
pbmc <- FindClusters(pbmc,resolution = 0.6)
DimPlot(pbmc,reduction = "umap",group.by = "RNA_snn_res.0.6",raster = F)
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
