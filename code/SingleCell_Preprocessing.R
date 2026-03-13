library(Seurat)
library(harmony)
library(ggplot2)
library(scop)
library(DoubletFinder)
reticulate::use_condaenv("F:/software/anaconda/envs/scop_env/")
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
  pbmc
})

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

seu_list<-lapply(seu_list,function(pbmc){
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
  pbmc<- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  scale.genes <-  VariableFeatures(pbmc)
  pbmc <- ScaleData(pbmc, features = scale.genes)
  pbmc<- RunPCA(pbmc, features = VariableFeatures(pbmc))
  pbmc<-RunUMAP(pbmc, dims = 1:20,reduction = 'pca')
  pbmc <- FindNeighbors(pbmc, dims = 1:20, reduction = "pca")
  pbmc <- FindClusters(pbmc,resolution = 0.6)
  sweep.res.list <- paramSweep(JoinLayers(pbmc), PCs = 1:20, sct = FALSE,num.cores = 16)
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
saveRDS(seu_list,"data/seu_list1.rds")

seu_list<-readRDS("data/seu_list1.rds")

pbmc<-merge(seu_list[[1]],seu_list[-1])

Idents(pbmc)<-"orig.ident"
#VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#pbmc<-scRNA
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc<- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
scale.genes <-  VariableFeatures(pbmc)
pbmc <- ScaleData(pbmc, features = scale.genes)
pbmc<- RunPCA(pbmc, features = VariableFeatures(pbmc))

seurat_merge_v5 <- IntegrateLayers(
  object = pbmc, 
  method = CCAIntegration, 
  orig.reduction = "pca", 
  new.reduction = "cca"
)
saveRDS(seurat_merge_v5,"data/cca1.rds")

# RPCA
seurat_merge_v5 <- IntegrateLayers(
  object = seurat_merge_v5, 
  method = RPCAIntegration, 
  orig.reduction = "pca", 
  new.reduction = "rpca"
)
saveRDS(seurat_merge_v5,"data/Rpca1.rds")

pbmc<-readRDS("data/Rpca1.rds")
#去批次----
scRNA_harmony <- RunHarmony(pbmc, group.by.vars = "orig.ident")
DimPlot(scRNA_harmony, reduction = "harmony", group.by = "orig.ident",raster = F)

ElbowPlot(scRNA_harmony,reduction = 'harmony',ndims = 50)
scRNA_harmony <- RunUMAP(scRNA_harmony, dims = 1:30,reduction = 'harmony')
#scRNA_harmony<-RunTSNE(scRNA_harmony,dims = 1:30,reduction = "harmony")
pbmc <- FindNeighbors(scRNA_harmony, dims = 1:30, reduction = "harmony")
pbmc <- FindClusters(pbmc,resolution = seq(0.1,0.1,0.1))

markers<-list(ProNeutrophil=c("DEFA3","CAMP","LCN2"),
              HSCs=c("CYTL1"),
              PlassmaCells=c("JCHAIN","MZB1","CD79A"),
              BCells=c("CD79B","MS4A1"),
              `T and NK`=c("CD3D","CD3E","CD3G"),
              Neutrophils=c("FCGR3B","CSF3R","S100A9"),
              Basophils=c("HDC","MS4A2","CPA3"),
              Monocyte=c("CD14","CD68","S100A12"),
              pDCs=c("CLEC4C","IL3RA","LILRA4"),
              Erythrocytes=c("HBB","HBA1","HBA2"),
              Platelets=c("PPBP","PF4","GP9"))

do_DotPlot(pbmc,features = markers,group.by = "RNA_snn_res.0.8")

celltype<-data.frame(cluster=pbmc$RNA_snn_res.0.8,celltype=NA)
celltype$celltype[celltype$cluster%in%c(23)]<-"ProNeutrophil"
celltype$celltype[celltype$cluster%in%c(33)]<-"HSCs"
celltype$celltype[celltype$cluster%in%c(21)]<-"PlassmaCells"
celltype$celltype[celltype$cluster%in%c(12, 10, 27, 31, 36)]<-"BCells"
celltype$celltype[celltype$cluster%in%c(7,2,3,1,9,29,16)]<-"TCells"
celltype$celltype[celltype$cluster%in%c(0,13,15,22,26)]<-"NK"
celltype$celltype[celltype$cluster%in%c(4,8,19,30,32)]<-"Neutrophils"
celltype$celltype[celltype$cluster%in%c(18)]<-"Basophils"
celltype$celltype[celltype$cluster%in%c(5,14, 6,25, 11, 17, 24, 35)]<-"Monocyte"
celltype$celltype[celltype$cluster%in%c(28)]<-"pDCs"
celltype$celltype[celltype$cluster%in%c(34)]<-"Erythrocytes"
celltype$celltype[celltype$cluster%in%c(20)]<-"Platelets"
pbmc$celltype<-celltype$celltype
pbmc$celltype<-factor(pbmc$celltype,levels = names(markers))
CellDimPlot(
  srt = pbmc,
  group.by = c("celltype"),
  reduction = "UMAP",
  theme_use = "theme_blank"
)
ggsave("result/singlecell/celltype.pdf",height = 6,width = 7)

do_DotPlot(pbmc,features = markers,group.by = "celltype")
ggsave("result/singlecell/cellanno_dot.pdf",width = 12,height = 6)

per_data<-as.data.frame(table(pbmc$celltype,pbmc$orig.ident))
sample<-unique(pbmc$orig.ident)
percent<-data.frame()
for(i in 1:length(sample)){
  tmp<-per_data[per_data$Var2==sample[i],]
  tmp$percent<-prop.table(tmp$Freq)
  percent<-rbind(percent,tmp)
}
percent$group<-gsub("-.*","",percent$Var2)

library(rstatix)
library(dplyr)
library(ggpubr)

state_result <- percent %>%
  group_by(Var1) %>%
  wilcox_test(percent ~ group, p.adjust.method = "BH")

state_result$p <- signif(state_result$p, 3) 
state_result$p2 <- ifelse(state_result$p < 0.05,
                          ifelse(state_result$p < 0.01,
                                 ifelse(state_result$p < 0.001,
                                        ifelse(state_result$p < 0.0001,"****","***"),"**"),"*"), "")


ggplot(percent,aes(Var1,percent))+
  geom_boxplot(aes(fill = group))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  stat_pvalue_manual(state_result, x="Var1", y.position=0.6, label="p2", size=4.5, color="red")+
  scale_fill_brewer(palette = "Set3")+
  labs(x="Celltype",y="Cell proportion")
ggsave("result/singlecell/cellproportion.pdf",width = 7,height = 6)
source("code/util.R")

NK<-pbmc[,pbmc$celltype=="NK"]
NK_res<-cell_DEG_enrich(NK,top = 15)
pDCs<-pbmc[,pbmc$celltype=="pDCs"]
pDCs_res<-cell_DEG_enrich(pDCs)
Platelets<-pbmc[,pbmc$celltype=="Platelets"]
Platelets_res<-cell_DEG_enrich(Platelets)
psa<-readRDS("result/singlecell/cell_per.rds")+ggtitle("A")+
  theme(plot.title.position = "plot")
psb<-NK_res[[1]]+ggtitle("B")
pse<-NK_res[[2]]+ggtitle("E")+theme(legend.key.height = unit(0.5,"cm"))
psc<-pDCs_res[[1]]+ggtitle("C")
psf<-pDCs_res[[2]]+ggtitle("F")+theme(legend.key.height = unit(0.5,"cm"))
psd<-Platelets_res[[1]]+ggtitle("D")
psg<-Platelets_res[[2]]+ggtitle("G")+theme(legend.key.height = unit(0.5,"cm"))
tmp1<-grid.arrange(psb,psc,psd,nrow=3)
tmp2<-wrap_plots(pse,psf,psg,nrow=3,heights=c(3,1,2))
tmp3<-cowplot::plot_grid(tmp1,tmp2,nrow=1,rel_widths=c(2,3))
dev.off()
pdf("result/singlecell/figS_Cel_DEG.pdf",width = 14,height = 14)
cowplot::plot_grid(psa,tmp3,nrow=2,rel_heights = c(1,3))
dev.off()
