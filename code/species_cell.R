library(ggplot2)
library(dplyr)
dif_species<-read.csv("D:/master/GGMP7009/result/dif_species/dif_genus.csv",row.names = 1)
dif_species$fdr<-p.adjust(dif_species$pvalue,method = "BH")
dif_genus<-dif_species[dif_species$fdr<0.05,]
dif_genus$feature<-gsub("g__","",dif_genus$feature)

lefse<-dif_genus
lefse$lda<-ifelse(lefse$enrich_group=="healthy",
                                  -lefse$ef_lda,
                  lefse$ef_lda)
marker_table<-lefse
ggplot(marker_table,aes(x=lda,y=reorder(feature,lda),fill=enrich_group))+
  geom_col()+
  theme_test()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(size = 16,face = "bold"),
        legend.text = element_text(size = 14),
        axis.title = element_text(size=14,face = "bold"),
        axis.text = element_text(size=10),
        plot.margin = margin(t=5,b=15,l=40,r=20),
        legend.position = "right",
        plot.title = element_text(hjust = -0.14,vjust = 0),
        title = element_text(size = 16,face="bold"))+
  geom_text(data = marker_table[marker_table$enrich_group == "healthy",], aes(y = feature, x = 0.1, label = feature),
            hjust = 0, size = 3.5) +
  geom_text(data = marker_table[marker_table$enrich_group == "hypertension",], aes(y = feature, x = -0.1, label = feature),
            hjust = 1, size = 3.5)+expand_limits(x=c(-6,5))+
  labs(x="LDA score",y="")+
  scale_fill_manual(values=c("#F9D5B2","#9FC3E2"))


bps<-read.table("result/scRPS_res/BPS_AUC.txt",sep = "\t",header = T,row.names = 1)
bps_p<-read.table("result/scRPS_res/pvalue_AUC.txt",sep = "\t",header = T,row.names = 1)
for(i in 1:ncol(bps_p)){
  bps_p[,i]<-p.adjust(bps_p[,i],method = "BH")
}

dif_genus_bps<-bps[,intersect(colnames(bps),dif_genus$feature)]
dif_genus_bps_p<-bps_p[,intersect(colnames(bps_p),dif_genus$feature)]
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

ggplot(bps_long,aes(celltype,bacteria,color = BPS_auc))+
  geom_point(aes(size = -log10(FDR)))+
  scale_color_gradient(low = "grey95",high = "slateblue2")+
  geom_text(aes(label = pvalue1),guide=F,color="black")+
  theme_test()+
  theme(axis.text.x = element_text(angle = 45,hjust=1))
