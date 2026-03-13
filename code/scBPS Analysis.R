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

magma_res<-read.table("data/magma_result.txt",header = T)
dif_genus<-read.csv("result/dif_species/dif_genus_intersect.csv")
magma_res<-magma_res[magma_res$TRAIT%in%dif_genus$genus,]
magma_list<-sapply(magma_res$GENESET, function(x){
  tmp<-strsplit(x,",|:")
  Z<-tmp[[1]][seq(2,length(tmp[[1]]),2)]|>as.numeric()
  names(Z)<-tmp[[1]][seq(1,length(tmp[[1]]),2)]
  sort(Z,decreasing = T)
})

kk_list<-sapply(magma_list,FUN = function(x){
  genes<-names(x)[1:100]
  genes<-bitr(genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  kk<-enrichKEGG(genes$ENTREZID)
  kk
})
names(kk_list)<-1:11
kk_list<-lapply(kk_list,function(x){
  x@result[x@result$pvalue<0.05,]
})

kk_name<-sapply(kk_list,function(x){
  x$Description
})|>unlist()
kk_freq<-table(kk_name)|>sort()|>as.data.frame()
kk_freq<-kk_freq[kk_freq$Freq>2,]
pb<-ggplot(kk_freq,aes(Freq,reorder(kk_name,Freq)))+
  geom_col(fill="skyblue",color="black")+
  theme_classic()+
  labs(y="",title = "B")+
  theme(
        legend.title = element_text(size = 16,face = "bold"),
        legend.text = element_text(size = 14),
        axis.title = element_text(size=14,face = "bold"),
        axis.text = element_text(size=10),plot.title.position = "plot",
        legend.position = "right",
        title = element_text(size = 16,face="bold"))

top100<-sapply(magma_list,FUN = function(x){
  names(x)[1:100]
})|>as.character()
top100<-top100[top100!=""]
top100<-as.data.frame(table(top100))
set.seed(123)
top100<-top100[order(top100$Freq,decreasing = T),]
wd<-ggplot(top100[1:200,], aes(label = top100, size = Freq)) +
  geom_text_wordcloud(aes(color=Freq),show.legend = F,
                      shape = "square",
                      eccentricity = 0.7) +
  theme_minimal()+
  scale_size_area(max_size = 4)+
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(name = "Spectral",n=11)))+
  theme(
    legend.title = element_text(size = 16,face = "bold"),
    legend.text = element_text(size = 14),
    axis.title = element_text(size=14,face = "bold"),
    axis.text = element_text(size=10),plot.title.position = "plot",
    legend.position = "right",
    title = element_text(size = 16,face="bold"))+
  labs(title = "A")
pdf("result/scBPS_res/enrichment(Fig2).pdf",width = 11,height = 5)
grid.arrange(wd,pb,nrow=1,widths=c(4,6))
dev.off()

bps<-read.table("result/scBPS_res/BPS_AUC.txt",sep = "\t",header = T,row.names = 1)
bps_p<-read.table("result/scBPS_res/pvalue_AUC.txt",sep = "\t",header = T,row.names = 1)
for(i in 1:ncol(bps_p)){
  bps_p[,i]<-p.adjust(bps_p[,i],method = "BH")
}

dif_genus_bps<-bps[,intersect(colnames(bps),dif_genus$genus)]
dif_genus_bps_p<-bps_p[,intersect(colnames(bps_p),dif_genus$genus)]
dif_genus_bps_p1<-dif_genus_bps_p
for(i in 1:ncol(dif_genus_bps_p)){
  dif_genus_bps_p[,i]<-case_when(dif_genus_bps_p[,i]<0.001~"***",
                                 dif_genus_bps_p[,i]<0.01~"**",
                                 dif_genus_bps_p[,i]<0.05~"*",
                                 .default = "")
}
pheatmap::pheatmap(dif_genus_bps,display_numbers = dif_genus_bps_p)

bps_long<-reshape2::melt(as.matrix(dif_genus_bps))
bps_plong<-reshape2::melt(as.matrix(dif_genus_bps_p))
bps_plong1<-reshape2::melt(as.matrix(dif_genus_bps_p1))

colnames(bps_long)<-c("celltype","bacteria","BPS_auc")
bps_long$pvalue1<-bps_plong$value
bps_long$FDR<-bps_plong1$value
bps_long$BPS_auc_strength<-ifelse(bps_long$BPS_auc>=0.3,ifelse(bps_long$BPS_auc>=0.5,">=0.5",">=0.3&<0.5"),"<0.3")
p3e<-ggplot(bps_long,aes(celltype,bacteria))+
  geom_point(aes(fill = -log10(FDR),
                 colour = BPS_auc_strength,
                 size=BPS_auc),stroke = 1.5,
             shape = 21)+
  scale_color_manual(values = c(">=0.5"="red",">=0.3&<0.5"="gold","<0.3"="grey90"))+
  scale_fill_gradient(low = "grey95",high = "slateblue2")+
  geom_text(aes(label = pvalue1),guide=F,color="black")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45,hjust=1),
        axis.text.y = element_text(color = ifelse(colnames(dif_genus_bps)%in%
                                                    dif_genus$genus[dif_genus$enrich_group=="hypertension"],
                                                  "red4","blue4")))
ggsave("result/scBPS_res/scBPS_dotplot.pdf",height = 5,width = 7)

pbmc<-readRDS("data/singlecell/pbmc_cellannoed.rds")
p3a<-CellDimPlot(pbmc,group.by = "RNA_snn_res.0.8",theme_use = "theme_blank",label = T)
p3b<-CellDimPlot(pbmc,group.by = "celltype",theme_use = "theme_blank")
markers<-list(ProNeutrophil=c("DEFA3","CAMP","LCN2"),
              HSCs=c("CYTL1"),
              PlassmaCells=c("JCHAIN","MZB1","CD79A"),
              BCells=c("CD79B","MS4A1"),
              TCells=c("CD3D","CD3E","CD3G"),
              NK=c("CD160", "KLRC1","NCAM1"),
              Neutrophils=c("FCGR3B","CSF3R","S100A9"),
              Basophils=c("HDC","MS4A2","CPA3"),
              Monocyte=c("CD14","CD68","S100A12"),
              pDCs=c("CLEC4C","IL3RA","LILRA4"),
              Erythrocytes=c("HBB","HBA1","HBA2"),
              Platelets=c("PPBP","PF4","GP9"))

p3c<-do_DotPlot(pbmc,features = markers,group.by = "celltype")

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


p3d<-ggplot(percent,aes(Var1,percent))+
  geom_boxplot(aes(fill = group))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  stat_pvalue_manual(state_result, x="Var1", y.position=0.6, label="p2", size=4.5, color="red")+
  scale_fill_manual(values = c("#9FC3E2","#F9D5B2"))+
  labs(x="Celltype",y="Cell proportion")
saveRDS(p3d,"result/singlecell/cell_per.rds")
# Extended analysis of monocytes ----
Monocyte<-pbmc[,pbmc$celltype=="Monocyte"]
sel_sp<-c("Bacteroides","Bifidobacterium","Butyricimonas","Desulfovibrio")
score<-read.table("result/scBPS_res/norm_score.tsv",sep = "\t",header = T,row.names = 1)
score_Mon<-score[colnames(Monocyte),sel_sp]
score_Mon<-reshape2::melt(as.matrix(score_Mon))
score_Mon$group<-Monocyte$group[match(score_Mon$Var1,colnames(Monocyte))]
p3f1<-ggplot(score_Mon,aes(y=group,x=value,fill = group))+
  geom_violin()+
  geom_boxplot(width=0.25,outlier.color = NA)+
  facet_grid(Var2~.)+
  ggpubr::stat_compare_means(comparisons = list(c("Control","Hypertension")),
                             method = "t.test",tip.length = 0,label.y = 8.5,label = "p.signif"
  )+
  scale_fill_manual(values = c("#9FC3E2","#F9D5B2"))+
  ggthemes::theme_clean()+
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        legend.position = "right")+
  labs(y="",x="BPS")
Tcel<-pbmc[,pbmc$celltype=="TCells"]
score_t<-score[colnames(Tcel),"Oscillospira",drop=F]
score_t<-reshape2::melt(as.matrix(score_t))
score_t$group<-Tcel$group[match(score_t$Var1,colnames(Tcel))]
p3f2<-ggplot(score_t,aes(y=group,x=value,fill = group))+
  geom_violin()+
  geom_boxplot(width=0.25,outlier.color = NA)+
  facet_grid(Var2~.)+
  ggpubr::stat_compare_means(comparisons = list(c("Control","Hypertension")),
                             method = "t.test",tip.length = 0,label.y = 6.5,label = "p.signif"
  )+
  scale_fill_manual(values = c("#9FC3E2","#F9D5B2"))+
  ggthemes::theme_clean()+
  theme(plot.background = element_blank(),
        legend.position = "none",
        legend.background = element_blank())+
  labs(y="",x="BPS")
p3f<-wrap_plots(p3f1,p3f2,heights = c(8,2))
mythemes<-theme(plot.title.position = "plot",
                plot.title = element_text(size = 14,face = "bold"))
p3a1<-p3a+ggtitle("A")+
  mythemes
p3b_1<-p3b+ggtitle("B")+
  mythemes
p3c_1<-p3c+ggtitle("C")+
  mythemes+
  guides(
    fill  = guide_colorbar(
      title.position = "top",direction = "vertical"
    )
  )+
  theme(legend.position = "right")
p3d1<-p3d+ggtitle("D")+
  mythemes
p3e1<-p3e+ggtitle("D")+
  mythemes
  #theme(
  #  legend.position = "bottom",
  #  legend.box = "horizontal"
  #)# +
  #guides(
  #  color = guide_legend(ncol = 1, title.position = "top"),
  #  size  = guide_legend(ncol = 1, title.position = "top"),
  #  fill  = guide_colorbar(
  #    title.position = "top",direction = "vertical"
  #  )
  #)
p3f1<-p3f
p3f1[[1]]<-p3f1[[1]]+ggtitle("E")+
  mythemes
tmp1<-grid.arrange(p3a1,p3b_1,nrow=1)
tmp2<-grid.arrange(tmp1,p3c_1,nrow=2,heights=c(2,2))
tmp3<-plot_grid(p3e1,p3f1,nrow=1,rel_widths = c(6,4))
#tmp4<-grid.arrange(tmp3,p3e1,ncol=2)
dev.off()
pdf("result/singlecell/Fig3.pdf",height = 14,width = 13)
grid.arrange(tmp2,tmp3,nrow=2,heights=c(6,4))
dev.off()

p3d1<-p3d+ggtitle("A")+
  mythemes

ps3b<-do_BarPlot(pbmc, 
           group.by = "celltype",
           split.by = "group",
           position = "fill",legend.position = "right",
           flip = F)+ggtitle("B")+
  mythemes
pdf("result/singlecell/Figs_per.pdf",width = 9,height = 5)
grid.arrange(p3d1,ps3b,nrow=1,widths=c(5,2.5))
dev.off()

BP<-read.csv("data/SraRunTable.csv",header = T)
BP$DBP<-gsub(".*/","",BP$Office.Blood.Pressure.mmHg.)|>as.numeric()
BP$SBP<-gsub("/.*","",BP$Office.Blood.Pressure.mmHg.)|>as.numeric()
Monocyte$SBP<-BP$SBP[match(Monocyte$orig.ident,BP$Library.Name)]
Monocyte$DBP<-BP$DBP[match(Monocyte$orig.ident,BP$Library.Name)]
Monocyte$Bacteroides_BPS<-score$Bacteroides[match(colnames(Monocyte),rownames(score))]
Monocyte$Bifidobacterium<-score$Bifidobacterium[match(colnames(Monocyte),rownames(score))]
Monocyte$Butyricimonas<-score$Butyricimonas[match(colnames(Monocyte),rownames(score))]
Monocyte$Desulfovibrio<-score$Desulfovibrio[match(colnames(Monocyte),rownames(score))]
Monocyte$SBP<-factor(Monocyte$SBP,levels = sort(unique(Monocyte$SBP)))

magma_res<-magma_res[magma_res$TRAIT%in%sel_sp,]
magma_list<-sapply(magma_res$GENESET, function(x){
  tmp<-strsplit(x,",|:")
  Z<-tmp[[1]][seq(2,length(tmp[[1]]),2)]|>as.numeric()
  names(Z)<-tmp[[1]][seq(1,length(tmp[[1]]),2)]
  sort(Z,decreasing = T)
})
top100<-sapply(magma_list,FUN = function(x){
  names(x)[1:100]
})|>as.character()
top100<-top100[top100!=""]
top100<-as.data.frame(table(top100))
wd<-ggplot(top100, aes(label = top100, size = Freq)) +
  geom_text_wordcloud(aes(color=Freq),show.legend = F) +
  theme_minimal()+
  scale_color_gradientn(colours = rev(c("#9E0142", "#D53E4F", "#F46D43" ,
                                        "#FDAE61" ,"gold", "#ABDDA4",
                                        "#66C2A5", "#3288BD", "#5E4FA2")))
wd



Monocyte<- FindVariableFeatures(Monocyte, selection.method = "vst", nfeatures = 2000)
scale.genes <-  VariableFeatures(Monocyte)
Monocyte <- ScaleData(Monocyte, features = scale.genes)
Monocyte<- RunPCA(Monocyte, features = VariableFeatures(Monocyte))
Monocyte <- RunHarmony(Monocyte, group.by.vars = "orig.ident")
CellDimPlot(Monocyte,reduction = "harmony",
            group.by = "orig.ident",
            theme_use = "theme_blank")
Monocyte <- RunUMAP(Monocyte, dims = 1:30,reduction = 'harmony')
Monocyte <- FindNeighbors(Monocyte, dims = 1:30, reduction = "harmony")
Monocyte <- FindClusters(Monocyte,resolution = 0.8)

p4a<-CellDimPlot(Monocyte,reduction = "umap",group.by = "seurat_clusters")
p4c<-FeatureDimPlot(Monocyte,features = c("FCGR3A","CD14","CST3","S100A8"))
p4b<-DotPlot(Monocyte,features = c("FCGR3A","CD14","CST3","S100A8"))

subcelltype<-data.frame(cluster=Monocyte$seurat_clusters,
                        subcelltype="Classical monocytes")
subcelltype$subcelltype[subcelltype$cluster%in%c(2,12)]<-"Non-classical monocytes"
subcelltype$subcelltype[subcelltype$cluster%in%c(5)]<-"Intermediate monocytes"
subcelltype$subcelltype[subcelltype$cluster%in%c(7)]<-"Non-classical–like"
Monocyte$subcelltype<-subcelltype$subcelltype
p4d<-CellDimPlot(Monocyte,reduction = "umap",group.by = "subcelltype")

BPS_df<-data.frame(BPS=c(Monocyte$Bacteroides_BPS,
                         Monocyte$Bifidobacterium,
                         Monocyte$Butyricimonas,
                         Monocyte$Desulfovibrio),
                   genus=rep(c("Bacteroides","Bifidobacterium","Butyricimonas","Desulfovibrio"),each=ncol(Monocyte)),
                   group=Monocyte$group,
                   subcelltype=subcelltype$subcelltype)
p4e<-ggplot(BPS_df,aes(subcelltype,BPS,fill=subcelltype))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.color = NA)+
  #geom_hline(yintercept = 1,linetype=2,color="red")+
  facet_grid(.~genus)+
  ggthemes::theme_clean()+
  scale_fill_brewer(palette = "Set2")+
  theme(plot.background = element_blank(),
        legend.background = element_blank(),
        axis.text.x = element_blank())
p4f<-ggplot(BPS_df,aes(group,BPS,fill=group))+
  geom_violin()+
  geom_boxplot(width=0.3,outlier.color = NA)+
  facet_grid(subcelltype~genus)+
  stat_compare_means(comparisons = list(c("Control","Hypertension")),
                     method = "t.test",label.y = 9,tip.length = 0)+
  ggthemes::theme_clean()+
  scale_fill_manual(values = c("#9FC3E2","#F9D5B2"))+
  theme(plot.background = element_blank(),
        legend.background = element_blank())

genesets<-getGmt("data/genesets/h.all.v2025.1.Hs.symbols.gmt")

cells_rankings <- AUCell_buildRankings(Monocyte[,Monocyte$subcelltype=="Classical monocytes"]@assays$RNA$data,
                                       nCores=8, plotStats=TRUE)

cells_AUC <- AUCell_calcAUC(genesets, cells_rankings, aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))

aucMatrix <- getAUC(cells_AUC)|>t()|>as.data.frame()


genesets_kegg<-getGmt("data/genesets/c2.all.v2025.1.Hs.symbols.gmt")
genesets_kegg<-genesets_kegg[grep("KEGG_",names(genesets_kegg))]

cells_rankings <- AUCell_buildRankings(Monocyte[,Monocyte$subcelltype=="Classical monocytes"]@assays$RNA$data,
                                       nCores=8, plotStats=TRUE)

cells_AUC <- AUCell_calcAUC(genesets_kegg, cells_rankings, aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))

kegg_aucMatrix <- getAUC(cells_AUC)|>t()|>as.data.frame()
M0<-Monocyte[,Monocyte$subcelltype=="Classical monocytes"]
cor_res<-list()
for(i in 1:4){
  tmp<-c()
  for(j in 1:ncol(kegg_aucMatrix)){
    cor<-cor.test(M0@meta.data[,c("Bacteroides_BPS","Bifidobacterium",
                                  "Butyricimonas","Desulfovibrio")[i]],kegg_aucMatrix[,j],method = "spearman")
    tmp<-c(tmp,cor$estimate)
  }
  cor_res[[i]]<-tmp
}


cor1<-pathway_cor(M0,sp = "Bacteroides_BPS",cor_res = cor_res[[1]],top=5,kegg_aucMatrix = kegg_aucMatrix)
cor2<-pathway_cor(M0,sp = "Bifidobacterium",cor_res = cor_res[[2]],top=5,kegg_aucMatrix = kegg_aucMatrix)
cor3<-pathway_cor(M0,sp = "Butyricimonas",cor_res = cor_res[[3]],top=10,kegg_aucMatrix = kegg_aucMatrix)

cor1_data<-cor1[[2]]
cor2_data<-cor2[[2]]
cor3_data<-cor3[[2]]

cor_all<-rbind(cor1_data,cor2_data)
cor_all$Var2<-as.character(cor_all$Var2)
cor1_data$Var2<-as.character(cor1_data$Var2)
cor_all$Var2[cor_all$Var2=="KEGG_MEDICUS_PATHOGEN_HCMV_UL33_TO_GNAQ_PLCB_G_CALCINEURIN_SIGNALING_PATHWAY"]<-
  "KEGG_MEDICUS_PATHOGEN_HCMV_UL33_TO_GNAQ\nPLCB_G_CALCINEURIN_SIGNALING_PATHWAY"
cor1_data$Var2[cor1_data$Var2=="KEGG_MEDICUS_PATHOGEN_HCMV_UL33_TO_GNAQ_PLCB_G_CALCINEURIN_SIGNALING_PATHWAY"]<-
  "KEGG_MEDICUS_PATHOGEN_HCMV_UL33_TO_GNAQ\nPLCB_G_CALCINEURIN_SIGNALING_PATHWAY"
p4g<-ggplot(cor_all,aes(x=Var1))+
  geom_line(data = cor1_data,mapping = aes(y=value,color=Var2),linewidth = 1.5)+
  geom_line(data = cor2_data,mapping = aes(y=value,color=Var2),linewidth = 1.5)+
  geom_point(aes(y=value,color=Var2,shape=genus),size=4)+
  labs(x="BPS (quintile bins)",
       y="Mean score of KEGG pathways")+
  theme_classic()+
  scale_color_brewer(name = "",palette = "Set3")+
  scale_x_continuous(breaks = seq(0,10,2))

cor3_data$Var2<-as.character(cor3_data$Var2)
cor3_data$Var2[cor3_data$Var2=="KEGG_MEDICUS_REFERENCE_ANTIGEN_PROCESSING_AND_PRESENTATION_BY_MHC_CLASS_II_MOLECULES"]<-
  "KEGG_MEDICUS_REFERENCE_ANTIGEN_PROCESSING\nAND_PRESENTATION_BY_MHC_CLASS_II_MOLECULES"
cor3_data$Var2[cor3_data$Var2=="KEGG_MEDICUS_PATHOGEN_HTLV_1_TAX_TO_NFY_MEDIATED_TRANSCRIPTION"]<-
  "KEGG_MEDICUS_PATHOGEN_HTLV_1_TAX_TO_NFY\nMEDIATED_TRANSCRIPTION"

magma_list<-sapply(magma_res[magma_res$TRAIT%in%c("Bacteroides","Bifidobacterium"),2],function(x){
  strsplit(x,",")
})
names(magma_list)<-c("Bacteroides","Bifidobacterium")
magma_list<-sapply(magma_list,function(x){
  gsub(":.*","",x)
})
ggvenn::ggvenn(magma_list)

p4a1<-p4a+ggtitle("A")+
  mythemes
p4b1<-p4b+ggtitle("B")+
  mythemes
p4c[[1]]<-p4c[[1]]+ggtitle("C")+
  mythemes
p4c
p4d1<-p4d+ggtitle("D")+
  mythemes
p4e1<-p4e+ggtitle("E")+
  mythemes
p4f1<-p4f+ggtitle("F")+
  mythemes+
  theme(axis.text.x = element_blank())
p4g1<-p4g+ggtitle("G")+
  mythemes+
  theme(legend.background = element_blank())

tmp1<-grid.arrange(p4a1,p4b1,nrow=1)
tmp2<-plot_grid(p4c,p4d1,nrow=1)
tmp3<-grid.arrange(tmp1,tmp2,nrow=2)
tmp4<-grid.arrange(tmp3,p4f1,nrow=1,widths=c(4,2.3))
tmp5<-grid.arrange(p4e1,p4g1,nrow=1,widths=c(4,6))
dev.off()
pdf("result/singlecell/Fig4.pdf",height = 11,width = 17)
grid.arrange(tmp4,tmp5,nrow=2,heights=c(5,1.8))
dev.off()
p5a<-ggplot(cor3_data,aes(x=Var1))+
  geom_line(data = cor3_data,mapping = aes(y=value,color=Var2),linewidth = 1.5)+
  geom_point(aes(y=value,color=Var2),size=4)+
  labs(x="BPS (quintile bins)",
       y="Mean score of KEGG pathways")+
  theme_classic()+
  scale_color_brewer(name = "",palette = "Set3")+
  scale_x_continuous(breaks = seq(0,10,2))


M1<-Monocyte[,Monocyte$subcelltype=="Intermediate monocytes"]
cells_rankings <- AUCell_buildRankings(M1@assays$RNA$data,
                                       nCores=8, plotStats=TRUE)

cells_AUC <- AUCell_calcAUC(genesets_kegg, cells_rankings, aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))

kegg_aucMatrix <- getAUC(cells_AUC)|>t()|>as.data.frame()
cor_res1<-c()
for(j in 1:ncol(kegg_aucMatrix)){
  cor<-cor.test(M1@meta.data[,c("Butyricimonas")],kegg_aucMatrix[,j],method = "spearman")
  cor_res1<-c(cor_res1,cor$estimate)
}

but_cor<-pathway_cor(M1,sp="Butyricimonas",cor_res1,top = 10,kegg_aucMatrix=kegg_aucMatrix)
but_cor[[1]]
but_cor1<-but_cor[[2]]
but_cor1$Var2<-as.character(but_cor1$Var2)
but_cor1$Var2[but_cor1$Var2=="KEGG_MEDICUS_PATHOGEN_HTLV_1_TAX_TO_NFY_MEDIATED_TRANSCRIPTION"]<-
  "KEGG_MEDICUS_PATHOGEN_HTLV_1_TAX_TO_NFY\nMEDIATED_TRANSCRIPTION"
but_cor1$Var2[but_cor1$Var2=="KEGG_MEDICUS_REFERENCE_ANTIGEN_PROCESSING_AND_PRESENTATION_BY_MHC_CLASS_II_MOLECULES"]<-
  "KEGG_MEDICUS_REFERENCE_ANTIGEN_PROCESSING_AND\nPRESENTATION_BY_MHC_CLASS_II_MOLECULES"
p5b<-ggplot(but_cor1,aes(x=Var1))+
  geom_line(data = but_cor1,mapping = aes(y=value,color=Var2),linewidth = 1.5)+
  geom_point(aes(y=value,color=Var2),size=4)+
  labs(x="BPS (quintile bins)",
       y="Mean score of KEGG pathways")+
  theme_classic()+
  scale_color_brewer(name = "",palette = "Set3")+
  scale_x_continuous(breaks = seq(0,10,2))


M2<-Monocyte[,Monocyte$subcelltype%in%c("Non-classical monocytes","Non-classical–like")]
cells_rankings <- AUCell_buildRankings(M2@assays$RNA$data,
                                       nCores=8, plotStats=TRUE)

cells_AUC <- AUCell_calcAUC(genesets_kegg, cells_rankings, aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))

kegg_aucMatrix <- getAUC(cells_AUC)|>t()|>as.data.frame()
cor_res2<-c()
for(j in 1:ncol(kegg_aucMatrix)){
  cor<-cor.test(M2@meta.data[,c("Butyricimonas")],kegg_aucMatrix[,j],method = "spearman")
  cor_res2<-c(cor_res2,cor$estimate)
}

but_cor2<-pathway_cor(M2,sp="Butyricimonas",cor_res2,top = 10,kegg_aucMatrix=kegg_aucMatrix )
but_cor2[[1]]
but_cor2<-but_cor2[[2]]
but_cor2$Var2<-as.character(but_cor2$Var2)
but_cor2$Var2[but_cor2$Var2=="KEGG_MEDICUS_VARIANT_MUTATION_INACTIVATED_PINK1_TO_ELECTRON_TRANSFER_IN_COMPLEX_I"]<-
  "KEGG_MEDICUS_VARIANT_MUTATION_INACTIVATED_PINK1\nTO_ELECTRON_TRANSFER_IN_COMPLEX_I"
but_cor2$Var2[but_cor2$Var2=="KEGG_MEDICUS_REFERENCE_MITOCHONDRIAL_COMPLEX_UCP1_IN_THERMOGENESIS"]<-
  "KEGG_MEDICUS_REFERENCE_MITOCHONDRIAL\nCOMPLEX_UCP1_IN_THERMOGENESIS"
but_cor2$Var2[but_cor2$Var2=="KEGG_MEDICUS_REFERENCE_ELECTRON_TRANSFER_IN_COMPLEX_I"]<-
  "KEGG_MEDICUS_REFERENCE_ELECTRON\nTRANSFER_IN_COMPLEX_I"
p5c<-ggplot(but_cor2,aes(x=Var1))+
  geom_line(data = but_cor2,mapping = aes(y=value,color=Var2),linewidth = 1.5)+
  geom_point(aes(y=value,color=Var2),size=4)+
  labs(x="BPS (quintile bins)",
       y="Mean score of KEGG pathways")+
  theme_classic()+
  scale_color_brewer(name = "",palette = "Set3")+
  scale_x_continuous(breaks = seq(0,10,2))



sub_Mono<-Monocyte[,Monocyte$group=="Control"&Monocyte$subcelltype%in%c("Classical monocytes",
                                                                        "Intermediate monocytes")]
sub_Mono <- RunSlingshot(
  srt = sub_Mono,
  group.by = "subcelltype",
  start="Classical monocytes",
  reduction = "umap"
)

sub_Monoh<-Monocyte[,Monocyte$group=="Hypertension"&Monocyte$subcelltype%in%c("Classical monocytes",
                                                                              "Intermediate monocytes")]
sub_Monoh <- RunSlingshot(
  srt = sub_Monoh,
  group.by = "subcelltype",
  start="Classical monocytes",
  reduction = "umap"
)

trace(CellDimPlot,edit = T) #526加入legend_list<-NULL
p5d<-CellDimPlot(sub_Monoh, group.by = "subcelltype", reduction = "umap", 
                 lineages = "Lineage1")
p5e<-FeatureDimPlot(
  sub_Monoh,
  features = paste0("Lineage", 1:3),
  reduction = "UMAP",
  theme_use = "theme_blank",
  lineages = "Lineage1"
)


cells_rankings <- AUCell_buildRankings(sub_Mono@assays$RNA$data,
                                       nCores=8, plotStats=TRUE)
ant_pro<-genesets_kegg[grep("ANTIGEN_PROCESSING",names(genesets_kegg))]
cells_AUC <- AUCell_calcAUC(ant_pro, cells_rankings, aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))

ant_pro_aucMatrix <- getAUC(cells_AUC)|>t()|>as.data.frame()
pseudotime<-sub_Mono@tools[["Slingshot_subcelltype_umap"]]@assays@data@listData[["pseudotime"]]
ant_pro_aucMatrix$pseudotime<-pseudotime[,1]|>as.numeric()
ggplot(ant_pro_aucMatrix,aes(pseudotime,KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION))+
  geom_point()


cells_rankings <- AUCell_buildRankings(sub_Monoh@assays$RNA$data,
                                       nCores=8, plotStats=TRUE)
cells_AUC <- AUCell_calcAUC(ant_pro, cells_rankings, aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))

ant_pro_aucMatrixh <- getAUC(cells_AUC)|>t()|>as.data.frame()
pseudotime<-sub_Monoh@tools[["Slingshot_subcelltype_umap"]]@assays@data@listData[["pseudotime"]]
ant_pro_aucMatrixh$pseudotime<-pseudotime[,1]|>as.numeric()

ant_pro_aucMatrix$group<-"Control"
ant_pro_aucMatrixh$group<-"Hypertension"
ant_pro_auc<-rbind(ant_pro_aucMatrix,ant_pro_aucMatrixh)


p5f<-ggplot(ant_pro_auc,aes(pseudotime,KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION))+
  #geom_point(aes(color=group))+
  geom_smooth(aes(color=group))+
  theme_classic()+
  scale_color_manual(values = c("blue4","red4"))+
  labs(y="KEGG_ANTIGEN_PROCESSING\nAND_PRESENTATION")

p5g<-ggplot(ant_pro_auc,aes(pseudotime,KEGG_MEDICUS_REFERENCE_ANTIGEN_PROCESSING_AND_PRESENTATION_BY_MHC_CLASS_II_MOLECULES))+
  #geom_point(aes(color=group))+
  geom_smooth(aes(color=group))+
  theme_classic()+
  scale_color_manual(values = c("blue4","red4"))+
  labs(y="KEGG_MEDICUS_REFERENCE_ANTIGEN\nPROCESSING_AND_PRESENTATION\nBY_MHC_CLASS_II_MOLECULES")


BP<-read.csv("data/SraRunTable.csv",header = T)
df<-ant_pro_aucMatrixh
df$sample<-sub_Monoh$orig.ident
df$BMI<-BP$BMI[match(df$sample,BP$Library.Name)]
df$Age<-BP$AGE[match(df$sample,BP$Library.Name)]
df$sex<-BP$sex[match(df$sample,BP$Library.Name)]
df$Lineage1<-sub_Monoh$Lineage1
df$Butyricimonas<-sub_Monoh$Butyricimonas
ml_g = lrn("regr.ranger", num.trees = 100, mtry = 2, min.node.size = 2, max.depth = 5)
ml_m = lrn("regr.ranger", num.trees = 100, mtry = 2, min.node.size = 2, max.depth = 5)
df$sex<-as.factor(df$sex)
df<-df[,c(10:12,14,13,1,5)]
df<-model.matrix(~.-1,df)
res<-data.frame()
for(i in 6:8){
  double_data<-df[,c(1:5,i)]
  double_data<-scale(double_data)|>as.data.frame()
  double_data<-na.omit(double_data)|>as.data.frame()
  model_data<-DoubleMLData$new(double_data,
                               y_col = colnames(df)[i],
                               d_cols = "Butyricimonas",
                               x_cols = c("Age","sexfemale","sexmale","BMI"))
  set.seed(123)
  model = DoubleMLPLR$new(model_data, ml_g, ml_m)
  model$fit()
  tmp<-data.frame(ATE=model$coef,lower=model$confint()[1],
                  upper=model$confint()[2],pvalue=model$pval)
  res<-rbind(res,tmp)
}
res$outcome<-colnames(df)[6:8]
res$outcome<-as.character(res$outcome)
res$outcome[res$outcome=="KEGG_MEDICUS_REFERENCE_ANTIGEN_PROCESSING_AND_PRESENTATION_BY_MHC_CLASS_II_MOLECULES"]<-
  "KEGG_MEDICUS_REFERENCE_ANTIGEN_PROCESSING\nAND_PRESENTATION_BY_MHC_CLASS_II_MOLECULES"
p5h<-double_vis(res,celltype = "Monocyte")

p5a1<-p5a+ggtitle("A")+
  mythemes
p5b1<-p5b+ggtitle("B")+
  mythemes
p5c1<-p5c+ggtitle("C")+
  mythemes
p5d1<-p5d+ggtitle("D")+
  mythemes
p5e1<-p5e+ggtitle("E")+
  mythemes
p5f1<-p5f+ggtitle("F")+
  mythemes
p5g1<-p5g+ggtitle("G")+
  mythemes
p5h1<-p5h+ggtitle("H")+
  mythemes

tmp1<-grid.arrange(p5a1,p5b1,nrow=1)
tmp2<-plot_grid(p5d1,p5e1,nrow=2)
tmp3<-grid.arrange(p5c1,tmp2,nrow=1,widths=c(5,3))
tmp4<-grid.arrange(tmp3,p5f1,nrow=1,widths=c(7,3))
tmp5<-grid.arrange(p5g1,p5h1,nrow=1,widths=c(3,6))
dev.off()
pdf("result/singlecell/Fig5.pdf",width = 16,height = 11)
grid.arrange(tmp1,tmp4,tmp5,nrow=3,heights=c(3,3,3))
dev.off()


pbmc$celltype1<-pbmc$celltype
pbmc$celltype<-as.character(pbmc$celltype)
Tcel<-readRDS("data/singlecell/Tcel_anno.rds")
pbmc$celltype[colnames(pbmc)%in%colnames(Tcel)[Tcel$subcelltype%in%c("CD4 Naive","CD4 TCM", "CD4 TEM","CD4 CTL","Treg")]]<-"CD4 T"
pbmc$celltype[colnames(pbmc)%in%colnames(Tcel)[Tcel$subcelltype%in%c("CD8 Naive","CD8 TCM", "CD8 TEM","CD8 CTL")]]<-"CD8 T"
pbmc$celltype[colnames(pbmc)%in%colnames(Tcel)[Tcel$subcelltype%in%c("dnT")]]<-"dnT"
pbmc$celltype[colnames(pbmc)%in%colnames(Tcel)[Tcel$subcelltype%in%c("MAIT")]]<-"MAIT"
pbmc$celltype[colnames(pbmc)%in%colnames(Monocyte)[Monocyte$subcelltype%in%c("Classical monocytes")]]<-"Classical monocytes"
pbmc$celltype[colnames(pbmc)%in%colnames(Monocyte)[Monocyte$subcelltype%in%c("Intermediate monocytes")]]<-"Intermediate monocytes"
pbmc$celltype[colnames(pbmc)%in%colnames(Monocyte)[Monocyte$subcelltype%in%c("Non-classical monocytes")]]<-"Non-classical monocytes"
pbmc$celltype[colnames(pbmc)%in%colnames(Monocyte)[Monocyte$subcelltype%in%c("Non-classical–like")]]<-"Non-classical–like monocytes"
cellchat.con<-cellchat_fun(pbmc)
cellchat.hy<-cellchat_fun(pbmc,group = "Hypertension")
object.list <- list(Control = cellchat.con, 
                    Hypertension = cellchat.hy) #对照组(NL)在前，比较组(LS)在后，注意顺序
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
p1<-compareInteractions(cellchat, show.legend = F, group = c(1,2))
p2<-compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p1<-p1+ggtitle("A")+
  mythemes+mythemes2
p2<-p2+
  mythemes+mythemes2
tmp1<-grid.arrange(p1,p2,nrow=2)

p1<-netVisual_diffInteraction(cellchat, weight.scale = T)
p2<-netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
p1<-as.ggplot(~replayPlot(p1))
p2<-as.ggplot(~replayPlot(p2))
p1<-p1+ggtitle("B")+mythemes
tmp2<-grid.arrange(p1,p2,nrow=1)
tmp2<-grid.arrange(tmp1,tmp2,nrow=1,widths=c(2,8))

p3 <- netVisual_heatmap(cellchat)
p4 <- netVisual_heatmap(cellchat, measure = "weight")
g <- grid.grabExpr(draw(p3))
p3 <- as.ggplot(g)
g <- grid.grabExpr(draw(p4))
p4 <- as.ggplot(g)
p3<-p3+ggtitle("C")+
  mythemes
tmp3<-grid.arrange(p3,p4,nrow=1)

pathways.show <- c("MHC-II") #选择目标信号通路
weight.max <- getMaxWeight(object.list, 
                           slot.name = c("netP"), 
                           attribute = pathways.show) #控制不同数据集的边权重
trace(netVisual_circle,edit = T)##111:igraph::E(g)$loop.angle<-0
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}


pathways.show <- c("MHC-II") 
par(mfrow = c(1,2), xpd = TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], 
                               signaling = pathways.show, 
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

g <- grid.grabExpr(draw(ht[[1]]))
ht1 <- as.ggplot(g)
g <- grid.grabExpr(draw(ht[[2]]))
ht2 <- as.ggplot(g)

ht1<-ht1+ggtitle("D")+
  mythemes
tmp4<-grid.arrange(ht1,ht2,nrow=T)
dev.off()

pdf("result/singlecell/figS_cellchat.pdf",width = 14,height = 16)
grid.arrange(tmp2,tmp3,tmp4,nrow=3)
dev.off()





hy_bar<-colnames(pbmc)[pbmc$group=="Hypertension"]
pbmc_hy<-pbmc[,pbmc$group=="Hypertension"]
pbmc_hy$Butyricimonas<-score[colnames(pbmc_hy),"Butyricimonas"]
pbmc_hy$Butyricimonas<-ifelse(pbmc_hy$Butyricimonas>median(pbmc_hy$Butyricimonas),"high BPS","low BPS")
cellchat.hy_h<-cellchat_fun(pbmc_hy[,pbmc_hy$Butyricimonas=="high BPS"],group = "Hypertension")
cellchat.hy_l<-cellchat_fun(pbmc_hy[,pbmc_hy$Butyricimonas=="low BPS"],group = "Hypertension")
object.list <- list(Low_BPS = cellchat.hy_l, 
                    High_BPS = cellchat.hy_h) #对照组(NL)在前，比较组(LS)在后，注意顺序
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
p1<-compareInteractions(cellchat, show.legend = F, group = c(1,2))
p2<-compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p1<-p1+ggtitle("A")+
  mythemes+mythemes2
p2<-p2+
  mythemes+mythemes2
tmp1<-grid.arrange(p1,p2,nrow=2)

p1<-netVisual_diffInteraction(cellchat, weight.scale = T)
p2<-netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
p1<-as.ggplot(~replayPlot(p1))
p2<-as.ggplot(~replayPlot(p2))
p1<-p1+ggtitle("B")+mythemes
tmp2<-grid.arrange(p1,p2,nrow=1)
tmp2<-grid.arrange(tmp1,tmp2,nrow=1,widths=c(2,8))

p3 <- netVisual_heatmap(cellchat)
p4 <- netVisual_heatmap(cellchat, measure = "weight")
g <- grid.grabExpr(draw(p3))
p3 <- as.ggplot(g)
g <- grid.grabExpr(draw(p4))
p4 <- as.ggplot(g)
p3<-p3+ggtitle("C")+
  mythemes
tmp3<-grid.arrange(p3,p4,nrow=1)

pathways.show <- c("MHC-II") #选择目标信号通路
weight.max <- getMaxWeight(object.list, 
                           slot.name = c("netP"), 
                           attribute = pathways.show) #控制不同数据集的边权重
trace(netVisual_circle,edit = T)##111:igraph::E(g)$loop.angle<-0
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "circle", 
                      edge.weight.max = weight.max[1], 
                      edge.width.max = 10, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}


pathways.show <- c("MHC-II") 
par(mfrow = c(1,2), xpd = TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], 
                               signaling = pathways.show, 
                               color.heatmap = "Reds",
                               title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

g <- grid.grabExpr(draw(ht[[1]]))
ht1 <- as.ggplot(g)
g <- grid.grabExpr(draw(ht[[2]]))
ht2 <- as.ggplot(g)

ht1<-ht1+ggtitle("D")+
  mythemes
tmp4<-grid.arrange(ht1,ht2,nrow=T)
dev.off()

pdf("result/singlecell/Fig7.pdf",width = 14,height = 16)
grid.arrange(tmp2,tmp3,tmp4,nrow=3)
dev.off()
