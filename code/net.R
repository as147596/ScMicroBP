library(ggplot2)
library(gghalves)
library(ggpubr)
library(NetCoMi)
library(igraph)
library(ggraph)
library(ggforce)
library(pheatmap)
library(tidyr)
otu<-read.table("data/GGMP7009/otu.tsv",row.names = 1,header = T,sep = "\t",comment.char = "",skip = 1)
taxa<-otu[,ncol(otu),drop=F]
meta<-read.csv("data/GGMP7009/meta.csv")
meta<-meta[!is.na(meta$DBP)&!is.na(meta$SBP),]
meta$group<-ifelse(meta$DBP<90&meta$SBP<140,"healthy","hypertension")
meta<-meta[which(meta$Antibiotics=="n"&meta$Medication=="n"),]
otu<-otu[,meta$SampleID]
taxa_tab<-sapply(taxa$taxonomy, function(x){
  strsplit(x,"; ")[[1]]
})|>unname()|>t()|>as.data.frame()
colnames(taxa_tab)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
taxa_tab$Species[taxa_tab$Species=="s__"]<-NA
taxa_tab$Genus[taxa_tab$Genus=="g__"]<-NA
taxa_tab$Family[taxa_tab$Family=="f__"]<-NA
taxa_tab$Order[taxa_tab$Order=="o__"]<-NA
taxa_tab$Class[taxa_tab$Class=="c__"]<-NA
taxa_tab$Phylum[taxa_tab$Phylum=="p__"]<-NA
rownames(taxa_tab)<-rownames(otu)
rownames(meta)<-meta$SampleID

otu$genus<-taxa_tab$Genus
genus<-aggregate(.~genus,otu,sum)
write.csv(genus,"data/GGMP7009/genus_df.csv",row.names=F,quote = F)

rm(otu)

rownames(genus)<-genus$genus
genus<-genus[,-1]
genus<-read.csv("data/GGMP7009/genus_df.csv",row.names = 1)
dif_genus<-read.csv("result/dif_species/dif_genus_intersect.csv",header=T)

keep<-apply(genus,1,function(x){
  mean(x!=0)>0.05
})
genus_df<-genus[keep,]

netobj<-netConstruct(data = t(genus_df),
                     dataType = "counts",
                     measure = "sparcc",
                     filtTax = "none",
                     filtSamp = "none",
                     normMethod = "TSS",
                     sparsMethod = "threshold",
                     thresh = 0,
                     seed = 123456)

sparcc_net<-netobj$assoMat1
diag(sparcc_net)<-0
res <- rm.get.threshold(sparcc_net, nr.thresholds = 50, 
                        interval = c(0.1, 0.3), discard.zeros = TRUE,
                        unfold.method = "spline")

netobj<-netConstruct(data = t(genus_df),
                     dataType = "counts",
                     measure = "sparcc",
                     filtTax = "none",
                     filtSamp = "none",
                     normMethod = "TSS",
                     sparsMethod = "threshold",
                     thresh = 0.15,
                     seed = 123456)

props <- netAnalyze(netobj, 
                    centrLCC = FALSE,
                    avDissIgnoreInf = TRUE,
                    sPathNorm = FALSE,
                    clustMethod = "cluster_walktrap",
                    hubPar = c("degree", "eigenvector"),
                    hubQuant = 0.9,
                    lnormFit = TRUE,
                    normDeg = FALSE,
                    normBetw = FALSE,
                    normClose = FALSE,
                    normEigen = FALSE)

plot(props, 
     nodeColor = "cluster", 
     nodeSize = "degree",
     repulsion = 0.8,
     rmSingles = TRUE,
     labelScale = FALSE,
     cexLabels = 0.2,
     nodeSizeSpread = 3,
     cexNodes = 2,
     hubBorderCol = "darkgray",
     title1 = "", 
     showTitle = TRUE,
     cexTitle = 2.3)
legend(0.35, 1.05, cex = 1.2, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)


net<-props$input$assoMat1
net<-reshape2::melt(net)
net<-net[net$value!=0,]
net$association<-ifelse(net$value>0,"+","-")
nodes<-data.frame(node=gsub("g__","",names(props$clustering$clust1)),
                  cluster=as.character(props$clustering$clust1))
dif_genus_control<-dif_genus[dif_genus[,2]=="healthy",]
dif_genus_hypertension<-dif_genus[dif_genus[,2]=="hypertension",]
nodes$enrich_group<-ifelse(nodes$node%in%dif_genus_control[,1],"healthy",ifelse(nodes$node%in%dif_genus_hypertension[,1],"hypertension","Notsignificant"))
net$Var1<-gsub("g__","",net$Var1)
net$Var2<-gsub("g__","",net$Var2)
p<-plot_graph(net,nodes)

saveRDS(p,"result/co_net/co_network.rds")
ggsave("result/co_net/co_network.pdf",plot = p,width = 9,height = 7)
per_df<-table(p$data$cluster,p$data$enrich_group)|>as.data.frame()
per_df<-pivot_wider(per_df,id_cols = Var1,names_from = Var2,values_from = Freq)|>as.data.frame()
per_df$size<-rowSums(per_df[,-1])
for(i in 2:ncol(per_df)){
  per_df[,i]<-paste0(per_df[,i],"(",
                     round(prop.table(per_df[,i])*100,2),
                     "%)")
}
colnames(per_df)<-c("cluster","enrich healty","enrich hypertension","Notsignificant","Size")
write.csv(per_df,"result/dif_species/cluster_percent.csv",row.names = F,quote = F)
degree_df<-sort(props[["centralities"]][["degree1"]],decreasing = T)
degree_df<-data.frame(nodes=gsub("g__","",names(degree_df)),
                      degree=unname(degree_df))[1:20,]
degree_df$type<-ifelse(degree_df$nodes%in%dif_genus[,1],"Significant","NotSignificant")
degree_df$type<-factor(degree_df$type,levels = c("Significant","NotSignificant"))
degree_df$nodes<-factor(degree_df$nodes,
                        levels = degree_df$nodes[order(degree_df$type,degree_df$degree)])
cluster<-props$clustering$clust1
c2<-gsub("g__","",names(cluster)[cluster==2])
#c4<-names(cluster)[cluster==4]
degree_df$cluster<-"other"
degree_df$cluster[degree_df$nodes%in%c2]<-"cluster2"
#degree_df$cluster[degree_df$nodes%in%c4]<-"cluster4"
p2<-ggplot(degree_df,aes(degree,y=nodes,fill=type))+
  geom_col()+
  geom_text(aes(x = degree+2,y=nodes,label = degree,color=cluster))+
  theme_test()+
  scale_fill_manual(name="",values = c("#F9D5B2","#9FC3E2"))+
  scale_color_manual(name="",values = c("red","grey"))+
  xlim(c(0,50))+
  theme(plot.background = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size=12))+
  scale_x_continuous(limits = c(0,40),breaks = seq(-0,40,10),labels = abs(seq(0,40,10)))+
  labs(x="",y="")
saveRDS(p2,"result/co_net/degree_top20.rds")
ggsave("result/co_net/degree_top20.pdf",width = 7,height = 5)
co_net_nodes<-nodes
c2_healthy<-co_net_nodes[co_net_nodes$cluster==2&co_net_nodes$enrich_group=="healthy",]
c2_hypertension<-co_net_nodes[co_net_nodes$cluster==2&co_net_nodes$enrich_group=="hypertension",]

genus_df_rel<-prop.table(as.matrix(genus_df),2)|>as.matrix()
rownames(genus_df_rel)<-gsub("g__","",rownames(genus_df_rel))
c2_health_df<-genus_df_rel[c2_healthy$node,]
c2_hy_df<-genus_df_rel[c2_hypertension$node,,drop=F]

box_data<-data.frame(value=c(colSums(c2_health_df),colSums(c2_hy_df)),
                     group=meta$group,
                     module=rep(c("c2_healthy","c2_hypertension"),each=ncol(genus_df)))
p3<-ggplot(box_data,aes(x=group,y=value))+
  geom_violin(aes(fill = group),alpha=0.7)+
  geom_boxplot(aes(fill = group),alpha=0.7,width=0.25)+
  facet_grid(~module)+
  ggpubr::stat_compare_means(
    comparisons = list(c("healthy","hypertension")),label.y = 0.9,
    method = "wilcox")+
  scale_fill_manual(values = c(hypertension="red4", healthy="#3C5488FF"))+
  theme_test()+
  theme(plot.background = element_blank())+
  labs(x="",y="sum of relative abundance")
ggsave("result/co_net/cluster2_sum_ab.pdf",height = 4,width = 7)
saveRDS(p3,"result/co_net/cluster2_sum_ab.rds")

