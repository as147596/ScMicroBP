library(Seurat)
library(harmony)
library(SCpubr)
library(ggplot2)
library(dplyr)
library(gghalves)
library(scop)
library(ggwordcloud)
library(clusterProfiler)
library(AUCell)
library(GSEABase)
library(org.Hs.eg.db)
library(monocle3)
library(SeuratWrappers)
library(slingshot)
library(RColorBrewer)
library(limma)
library(ggrepel)
library(ggpubr)
library(rstatix)
library(DoubleML)
library(mlr3)
library(mlr3learners)
source("code/util.R")
reticulate::use_condaenv("F:/software/anaconda/envs/scop_env/")
pbmc<-readRDS("data/singlecell/pbmc_cellannoed.rds")
Tcel<-pbmc[,pbmc$celltype=="TCells"]
rm(pbmc)

cell_markers <- list()

# CD4+ T细胞亚群
cell_markers$`Naive T cells` <- c("CCR7", "SELL", "TCF7", "IL7R", "GPR183", "LEF1", "LTB")
cell_markers$`Central Memory T cells` <- c("CCR7", "HSPA6", "MT1E", "MT1F", "CLDND1", "CD69", "GPR183", "IL7R", "KLF2", "TOB1")
cell_markers$`Effector Memory T cells` <- c("RPS19", "DUSP2", "CD44")
cell_markers$`Th1-like T cells` <- c("UCP2", "APCB1B", "TBX21", "IFNG")
cell_markers$`Th17-like T cells` <- c("LINC00513", "AC016831.7", "ADAM19", "XIST", "RORA", "IL17A", "IL17F", "RORC", "IL23R", "CTSH", "CAGC")
cell_markers$`Tfh` <- c("CXCL13", "CTLA4", "BATE", "PDCD1", "ICOS", "TNFRSF8", "BCL6", "TOX")
cell_markers$`Treg` <- c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "TNFRSF18", "MAGEH1", "SAT1", "CCR8", "KZF2", "IL10", "BATF")

# CD8+ T细胞亚群
cell_markers$`Naive T cells` <- c("CCR7", "SELL", "S100A8", "CST3", "AC020916.1", "IL7R", "LEF1", "TCF7")
cell_markers$`Central Memory T cells` <- c("CLND1", "RPS26", "GZMK", "CD44", "EOMES", "CD28", "CCR7", "DKK3")
cell_markers$`Effector Memory T cells` <- c("DUSP2", "GZMK", "CX3CR1")
cell_markers$`Teff` <- c("FGFBP2", "CX3CR1", "FCGR3A", "XLRG1", "IFNG", "GZMB", "GZMH", "PRF1", "NKG7", "GNLY")
cell_markers$`Trm` <- c("XCL1", "XCL2", "IL7R", "PRDM1", "TGFBR2", "ITGAL", "SELL")
cell_markers$`Tex` <- c("CXCL13", "TIGIT", "CTLA4", "PDCD1", "LAG3", "HAVCR2", "ENTPD1", "TOX", "TOX2", "LAYN", "TNFRSF9", "TCF7")
cell_markers$`Temma` <- c("FGFBP2", "ADGRG1", "FCGR3A")
cell_markers$`Cytotoxic` <- c("NKG7", "PRF1", "GZMK", "GZMA", "GZMB", "GZMH", "IFNG")
cell_markers$`TSTR` <- c("HSPA1A", "HSPA1B")
cell_markers$`TSEN` <- c("CD27")
cell_markers$`p-TEX` <- c("TCF7", "CD27", "CD28", "EOMES")

# 非常规T细胞亚群
cell_markers$`NKT cells` <- c("EOMES", "XCL1", "XCL2", "CXCR6", "TIGIT", "LAG3")
cell_markers$`MAIT cells` <- c("TRAV1-2", "SLC4A10", "GZMK", "KLRB1", "RORC", "RORA")
cell_markers$`IEL` <- c("CD160", "KIR2DL4", "TMIGD2", "ITGAE")
cell_markers$`DNT` <- c("GZMK")

Tcel<- FindVariableFeatures(Tcel, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(Tcel)
Tcel <- ScaleData(Tcel, features = scale.genes)
Tcel<- RunPCA(Tcel, features = VariableFeatures(Tcel))
Tcel <- RunHarmony(Tcel, group.by.vars = "orig.ident")
CellDimPlot(Tcel,reduction = "harmony",
            group.by = "orig.ident",
            theme_use = "theme_blank")
Tcel <- RunUMAP(Tcel, dims = 1:30,reduction = 'harmony')
Tcel <- FindNeighbors(Tcel, dims = 1:30, reduction = "harmony")
Tcel <- FindClusters(Tcel,resolution = 0.8)
p6a<-CellDimPlot(Tcel,reduction = "umap",group.by = "seurat_clusters",label = T)

cellmarker1<-list(CD4T=c("CD4", "CD3D","CD3E","CD3G"),
                  CD8T=c("CD8A", "CD8B"))


do_DotPlot(Tcel,features = cellmarker1)
FeatureDimPlot(Tcel,features = c("CD4","CD8A","CD"))
anno<-c("CD4T","CD8T","CD8T","CD4T","CD8T","CD4T","CD4T","CD4T","CD8T",
        "dnT","CD4T","dnT","CD4T","CD8T","CD4T","CD8T","CD4T")
names(anno)<-0:16
Tcel$subcelltype1<-anno[match(Tcel$seurat_clusters,names(anno))]|>unname()

CellDimPlot(Tcel,group.by = "subcelltype1")
cellmarker2<-list("CD4 Naive"=c("CD4","CCR7","IL7R","SELL","TCF7","LEF1"),
                  "CD4 TCM"=c("CCR7", "IL7R", "LTB"),
                  "CD4 TEM"=c("GZMK","S100A4"),
                  "CD4 CTL"=c("PRF1", "GZMB","NKG7"),
                  Treg=c("FOXP3", "IL2RA","CTLA4"),
                  "MAIT"=c("TRAV1-2","SLC4A10","KLRB1"),
                  "CD8 Naive"=c("CD8A", "CD8B","CCR7", "SELL", "IL7R", "TCF7","LEF1"),
                  "CD8 TCM"=c("CCR7", "IL7R", "LTB","CLND1", "GZMK", "CD44", "EOMES", "CD28", "CCR7", "DKK3"),
                  "CD8 TEM"=c("GZMK", "IFNG","DUSP2", "CX3CR1"),
                  "CD8 CTL"=c("GZMB", "PRF1", "NKG7"),
                  "CD8 Exhausted"=c("PDCD1", "LAG3", "HAVCR2")
)
do_DotPlot(Tcel,features = cellmarker2)
ggsave("result/singlecell/Tcell_dot.pdf",height = 10,width = 30)

subcelltype<-data.frame(cluster=Tcel$seurat_clusters,
                        subcelltype="CD4 Naive")
subcelltype$subcelltype[subcelltype$cluster%in%c(3,5)]<-"CD4 TCM"
subcelltype$subcelltype[subcelltype$cluster%in%c(4,6)]<-"CD4 TEM"
subcelltype$subcelltype[subcelltype$cluster%in%c(7)]<-"CD4 CTL"
subcelltype$subcelltype[subcelltype$cluster%in%c(10)]<-"Treg"
subcelltype$subcelltype[subcelltype$cluster%in%c(9)]<-"MAIT"
subcelltype$subcelltype[subcelltype$cluster%in%c(12)]<-"dnT"
subcelltype$subcelltype[subcelltype$cluster%in%c(2)]<-"CD8 Naive"
subcelltype$subcelltype[subcelltype$cluster%in%c(4)]<-"CD8 TEM"
subcelltype$subcelltype[subcelltype$cluster%in%c(5)]<-"CD8 TCM"
subcelltype$subcelltype[subcelltype$cluster%in%c(1,8,13,15)]<-"CD8 CTL"
Tcel$subcelltype<-factor(subcelltype$subcelltype,levels = c("CD4 Naive","CD4 TCM","CD4 TEM","CD4 CTL",
                                                            "Treg","MAIT","CD8 Naive","CD8 TCM",
                                                            "CD8 TEM","CD8 CTL","dnT"))
p6a<-CellDimPlot(Tcel,group.by = "seurat_clusters",label = T)
p6b<-CellDimPlot(Tcel,group.by = "subcelltype")
p6c<-do_DotPlot(Tcel,features = cellmarker2,group.by = "subcelltype")
saveRDS(Tcel,"data/singlecell/Tcel_anno.rds")

per_data<-as.data.frame(table(Tcel$subcelltype,Tcel$orig.ident))
sample<-unique(Tcel$orig.ident)
percent<-data.frame()
for(i in 1:length(sample)){
  tmp<-per_data[per_data$Var2==sample[i],]
  tmp$percent<-prop.table(tmp$Freq)
  percent<-rbind(percent,tmp)
}
percent$group<-gsub("-.*","",percent$Var2)
state_result <- percent %>%
  group_by(Var1) %>%
  wilcox_test(percent ~ group, p.adjust.method = "BH")

state_result$p <- signif(state_result$p, 3) 
state_result$p2 <- ifelse(state_result$p < 0.05,
                          ifelse(state_result$p < 0.01,
                                 ifelse(state_result$p < 0.001,
                                        ifelse(state_result$p < 0.0001,"****","***"),"**"),"*"), "")


p6d<-ggplot(percent,aes(Var1,percent))+
  geom_boxplot(aes(fill = group))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  stat_pvalue_manual(state_result, x="Var1", y.position=0.4, label="p2", size=4.5, color="red")+
  scale_fill_manual(values = c("#9FC3E2","#F9D5B2"))+
  labs(x="Celltype",y="Cell proportion")
score<-read.table("result/scBPS_res/norm_score.tsv",sep = "\t",header = T,row.names = 1)
Tcel$Oscillospira<-score$Oscillospira[match(colnames(Tcel),rownames(score))]
BPS_df<-data.frame(BPS=Tcel$Oscillospira,
                   genus="Oscillospira",
                   group=Tcel$group,
                   subcelltype=Tcel$subcelltype)
p6e<-ggplot(BPS_df,aes(subcelltype,BPS,fill=subcelltype))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.color = NA)+
  geom_hline(yintercept = 0.5,linetype=2,color="red")+
  facet_grid(.~genus)+
  ggthemes::theme_clean()+
  scale_fill_brewer(palette = "Set3")+
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45,hjust=1))

dnT<-Tcel[,Tcel$subcelltype=="dnT"]

score_tcel<-score[colnames(dnT),"Oscillospira",drop=F]
score_tcel$group<-dnT$group[match(rownames(score_tcel),colnames(dnT))]
p6f<-ggplot(score_tcel,aes(group,Oscillospira,fill = group))+
  geom_violin()+
  geom_boxplot(width=0.25,outlier.color = NA)+
  ggpubr::stat_compare_means(comparisons = list(c("Control","Hypertension")),
                             method = "t.test",tip.length = 0,label.y = 7
  )+
  scale_fill_manual(values = c("#9FC3E2","#F9D5B2"))+
  ggthemes::theme_clean()+
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45,hjust=1))+
  labs(x="",y="BPS")

genesets<-getGmt("data/genesets/h.all.v2025.1.Hs.symbols.gmt")

cells_rankings <- AUCell_buildRankings(Tcel@assays$RNA$data,
                                       nCores=8, plotStats=TRUE)

cells_AUC <- AUCell_calcAUC(genesets, cells_rankings, aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))

aucMatrix <- getAUC(cells_AUC)|>t()|>as.data.frame()
tmp<-c()
p<-c()
for(j in 1:ncol(aucMatrix)){
  cor<-cor.test(as.numeric(as.character(dnT@meta.data$Oscillospira)),aucMatrix[colnames(dnT),j],method = "spearman")
  tmp<-c(tmp,cor$estimate)
  p<-c(p,cor$p.value)
}
cor_res<-data.frame(geneset=colnames(aucMatrix),cor=tmp,p_value=p,genus="Oscillospira")
pathway_cor_res<-pathway_cor(dnT,sp = "Oscillospira",cor_res = tmp$cor,top=3,sel=c(
  "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_IL2_STAT5_SIGNALING"
),kegg_aucMatrix = aucMatrix[colnames(dnT),],ylab = "Mean score of HALLMARK pathway")
write.csv(cor_res,"result/singlecell/Tcell_pathway_cor.csv",quote = F,row.names = F)
p6g<-pathway_cor_res[[1]]
dnT_dep<-DEP(Tcel,"dnT",aucMatrix)
p6h<-dnT_dep[[2]]
cd4naive_dep<-DEP(Tcel,"CD4 Naive",aucMatrix)
cd4tcm_dep<-DEP(Tcel,"CD4 TCM",aucMatrix)
cd8naive_dep<-DEP(Tcel,"CD8 Naive",aucMatrix)
cd8tcm_dep<-DEP(Tcel,"CD8 TCM",aucMatrix)
Treg_dep<-DEP(Tcel,"Treg",aucMatrix)
p6i<-multi_vol(x1=cd4naive_dep[[1]],
          x2=cd4tcm_dep[[1]],
          x3=cd8naive_dep[[1]],
          x4=cd8tcm_dep[[1]],
          cells=c("CD4 Naive","CD4 TCM","CD8 Naive","CD8 TCM"))


CD4<-Tcel[,Tcel$subcelltype%in%c("CD4 Naive","Treg","dnT")&
            Tcel$group=="Control"]
CD4<-RunSlingshot(
  srt = CD4,
  group.by = "subcelltype",
  start="CD4 Naive",
  reduction = "UMAP"
)
FeatureDimPlot(
  CD4,
  features = paste0("Lineage", 1:3),
  reduction = "UMAP",
  theme_use = "theme_blank"
)


CD4h<-Tcel[,Tcel$subcelltype%in%c("CD4 Naive","Treg","dnT")&
            Tcel$group=="Hypertension"]
CD4h<-RunSlingshot(
  srt = CD4h,
  group.by = "subcelltype",
  start="CD4 Naive",
  reduction = "UMAP"
)
FeatureDimPlot(
  CD4h,
  features = paste0("Lineage", 1:3),
  reduction = "UMAP",
  theme_use = "theme_blank"
)


dnT$IL2_STAT5<-ifelse(aucMatrix[colnames(dnT),"HALLMARK_IL2_STAT5_SIGNALING"]>median(aucMatrix[colnames(dnT),"HALLMARK_IL2_STAT5_SIGNALING"]),
                "high","low")
dnT<-FindVariableFeatures(dnT, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(dnT)
dnT <- ScaleData(dnT, features = scale.genes)
dnT<- RunPCA(dnT, features = VariableFeatures(dnT))
dnT <- RunHarmony(dnT, group.by.vars = "orig.ident")
CellDimPlot(dnT,reduction = "harmony",
            group.by = "orig.ident",
            theme_use = "theme_blank")
dnT <- RunUMAP(dnT, dims = 1:30,reduction = 'harmony')
dnT <- FindNeighbors(dnT, dims = 1:30, reduction = "harmony")
dnT <- FindClusters(dnT,resolution = 0.8)
CellDimPlot(dnT,reduction = "umap",group.by = "seurat_clusters",label = T)
p6j<-CellDimPlot(dnT,reduction = "umap",group.by = "IL2_STAT5",palette ="Set1",theme_use = "theme_blank")
dnT<-RunSlingshot(
  srt = dnT,
  group.by = "IL2_STAT5",
  start="low",
  reduction = "UMAP",
  lineage_palette  ="Set1"
)
p6k<-FeatureDimPlot(
  dnT,
  features = paste0("Lineage", 1),
  reduction = "UMAP",
  theme_use = "theme_blank"
)

ggplot(dnT@meta.data,aes(Lineage1,Oscillospira))+
  geom_point(fill="slategray2",shape = 21,alpha=0.5)+
  geom_smooth(method = "lm",color="skyblue3",fill="skyblue3")+
  stat_cor(label.x = 1.5, label.y = 7.5, label.sep = "\n")+
  ggthemes::theme_base()+
  theme(plot.background = element_blank())

BP<-read.csv("data/SraRunTable.csv",header = T)
dnT$BMI<-BP$BMI[match(dnT$orig.ident,BP$Library.Name)]
dnT$IL2_STAT5<-aucMatrix[colnames(dnT),"HALLMARK_IL2_STAT5_SIGNALING"]
dnT$INFLAMMATORY_RESPONSE<-aucMatrix[colnames(dnT),"HALLMARK_INFLAMMATORY_RESPONSE"]
dnT$INTERFERON_ALPHA_RESPONSE<-aucMatrix[colnames(dnT),"HALLMARK_INTERFERON_ALPHA_RESPONSE"]
df<-dnT@meta.data[,c(4,5,76,72,73,74,77,78)]
ml_g = lrn("regr.ranger", num.trees = 100, mtry = 2, min.node.size = 2, max.depth = 5)
ml_m = lrn("regr.ranger", num.trees = 100, mtry = 2, min.node.size = 2, max.depth = 5)
df$sex<-as.factor(df$sex)
df<-model.matrix(~.-1,df)
res<-data.frame()
for(i in 6:9){
  double_data<-df[,c(1:5,i)]
  double_data<-scale(double_data)|>as.data.frame()
  double_data<-na.omit(double_data)|>as.data.frame()
  model_data<-DoubleMLData$new(double_data,
                               y_col = colnames(df)[i],
                               d_cols = "Oscillospira",
                               x_cols = c("age","sexfemale","sexmale","BMI"))
  set.seed(123)
  model = DoubleMLPLR$new(model_data, ml_g, ml_m)
  model$fit()
  tmp<-data.frame(ATE=model$coef,lower=model$confint()[1],
                  upper=model$confint()[2],pvalue=model$pval)
  res<-rbind(res,tmp)
}
res$outcome<-colnames(df)[6:9]

p6l<-double_vis(res)

p6a1<-p6a+ggtitle("A")+
  mythemes
p6b1<-p6b+ggtitle("B")+
  mythemes
p6c1<-p6c+ggtitle("C")+
  mythemes+
  theme(legend.position = "right",
        legend.background = element_blank(),
        legend.key.height = unit(0.5,"cm"))+
  guides(
    fill  = guide_colorbar(
      title.position = "top",direction = "vertical"
    )
  )
p6d1<-p6d+ggtitle("D")+
  mythemes
p6e1<-p6e+ggtitle("E")+
  mythemes
p6f1<-p6f+ggtitle("F")+
  mythemes
p6g1<-p6g+ggtitle("G")+
  mythemes+labs(y="Mean score of\nHallMARK pathway")
p6h1<-p6h+ggtitle("H")+
  mythemes
p6i1<-p6i+ggtitle("I")+
  mythemes
p6j1<-p6j+ggtitle("J")+
  mythemes
p6k1<-p6k+ggtitle("K")+
  mythemes
p6l1<-p6l+ggtitle("L")+
  mythemes

tmp1<-grid.arrange(p6a1,p6b1,p6d1,p6e1,nrow=1,widths=c(3,3,3,2))
tmp2<-grid.arrange(p6c1,p6f1,nrow=1,widths=c(7,1))
tmp3<-grid.arrange(p6g1,p6h1,nrow=2,heights=c(2,3))
tmp4<-grid.arrange(tmp3,p6i1,nrow=1,widths=c(2,3))
tmp5<-grid.arrange(p6j1,p6k1,p6l1,nrow=1,widths=c(3,3,7))
dev.off()
pdf("result/singlecell/Fig6.pdf",height = 16,width = 15)
grid.arrange(tmp1,tmp2,tmp4,tmp5,ncol=1,heights=c(2,2,3,1.5))
dev.off()
