library(MicrobiotaProcess)
library(microbiomeMarker)
library(ggplot2)
library(gghalves)
library(ggpubr)
library(curatedMetagenomicData)
library(MMUPHin)
library(UpSetR)
library(patchwork)
library(gridExtra)
library(grid)
source("code/util.R")
otu<-read.table("data/GGMP7009/otu.tsv",row.names = 1,header = T,sep = "\t",comment.char = "",skip = 1)
taxa<-otu[,ncol(otu),drop=F]
meta<-read.csv("data/GGMP7009/meta.csv")
meta<-meta[!is.na(meta$DBP)&!is.na(meta$SBP),]
meta$group<-ifelse(meta$DBP<90&meta$SBP<130,"healthy","hypertension")
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
mpse <- mp_import_dada2(seqtab=t(otu), taxatab=taxa_tab, sampleda=meta)
mpse %<>% mp_rrarefy()
mpse %<>% 
  mp_cal_alpha(.abundance=RareAbundance)
mpse %>% 
  mp_plot_alpha(
    .group=group, 
    .alpha=c(Observe, Chao1, ACE, Shannon, Simpson, Pielou)
  ) +
  scale_fill_manual(values=c("#00A087FF", "#3C5488FF"), guide="none") +
  scale_color_manual(values=c("#00A087FF", "#3C5488FF"), guide="none")
ggsave("result/dif_species/alpha.pdf",width = 6,height = 4)

mpse %<>%
  mp_cal_abundance( # for each samples
    .abundance = RareAbundance
  ) %>%
  mp_cal_abundance( # for each groups 
    .abundance=RareAbundance,
    .group=group
  )
mpse %<>% 
  mp_decostand(.abundance=Abundance)
mpse %<>% mp_cal_dist(.abundance=hellinger, distmethod="bray")
mpse %>% mp_plot_dist(.distmethod = bray, .group = group)
mpse %>% mp_plot_dist(.distmethod = bray, .group = group, group.test=TRUE, textsize=2)
mpse %<>% 
  mp_cal_pcoa(.abundance=hellinger, distmethod="bray")
mpse %<>%
  mp_adonis(.abundance=hellinger, .formula=~group, distmethod="bray", permutations=9999, action="add")
mpse %>% mp_extract_internal_attr(name=adonis)
mpse %>% 
  mp_plot_ord(
    .ord = pcoa, 
    .group = group, 
    .color = group, 
    .size = Observe, 
    .alpha = Shannon,
    ellipse = TRUE,
    show.legend = FALSE # don't display the legend of stat_ellipse 
  ) +
  scale_fill_manual(
    values = c("#00A087FF", "#3C5488FF"), 
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_color_manual(
    values=c("#00A087FF", "#3C5488FF"),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  ) +
  scale_size_continuous(
    range=c(0.5, 3),
    guide = guide_legend(keywidth=0.6, keyheight=0.6, label.theme=element_text(size=6.5))
  )
ggsave("result/dif_species/beat.pdf",width = 6,height = 6)
mpse %<>%
  mp_diff_analysis(
    .abundance = RelRareAbundanceBySample,
    .group = group,ldascore = 0,
    first.test.alpha = 0.05
  )
# The result is stored to the taxatree or otutree slot, you can use mp_extract_tree to extract the specific slot.
taxa.tree <- mpse %>% 
  mp_extract_tree(type="taxatree")
dif<-taxa.tree %>% select(label, nodeClass, LDAupper, LDAmean, LDAlower, Sign_group, pvalue, fdr) %>% dplyr::filter(!is.na(fdr))
write.csv(dif,"result/dif_species/difspecies.csv")


phy<-as.phyloseq(mpse)

mm_lefse <- run_lefse(
  phy,taxa_rank = "Genus",
  wilcoxon_cutoff = 0.05,
  group = "group",
  kw_cutoff = 0.05,
  multigrp_strat = T,
  lda_cutoff = 2
)

dif_genus_df<-mm_lefse@marker_table|>as.data.frame()|>as.matrix()

write.csv(dif_genus_df,"result/dif_species/dif_genus.csv")

dif_genus_pro<-read.csv("result/dif_species/dif_genus.csv",row.names = 1)
#dif_genus_test_df<-read.csv("data/test_dif_genus/genus_dif.csv")
dif_genus_pro$feature<-gsub("g__|_f__.*|_o__.*","",dif_genus_pro$feature)
dif_genus_p1<-draw_lefse(dif_genus_pro)

## metagenome ----

hypertension_meta<-sampleMetadata[which(sampleMetadata$disease=='hypertension'),]
healthy_meta<-hypertension_meta[which(hypertension_meta$age_category=='adult'),]
healthy_meta<-sampleMetadata[which(sampleMetadata$disease=='healthy'),]
healthy_meta<-healthy_meta[which(healthy_meta$age_category=='adult'),]
healthy_meta<-healthy_meta[which(healthy_meta$country=='CHN'|healthy_meta$country=='ITA'|healthy_meta$country=='AUT'),]
healthy_meta<-healthy_meta[which(healthy_meta$BMI>18.5&healthy_meta$BMI<23.9),]
healthy_meta<-healthy_meta[which(healthy_meta$age>55&healthy_meta$age<70),]
hypertension_meta<-rbind(hypertension_meta,healthy_meta)
ab_hypertension<-hypertension_meta |>
  dplyr::filter(body_site == 'stool') |>
  returnSamples("relative_abundance",rownames = "short")
abundance_hypertension<-ab_hypertension@assays@data@listData[["relative_abundance"]]
write.csv(abundance_hypertension,"data/abundance_short.csv")

ab_hypertension<-hypertension_meta |>
  dplyr::filter(body_site == 'stool') |>
  returnSamples("relative_abundance",rownames = "long")
#saveRDS(ab_hypertension,"data/rel_ab.rds")

ab_hypertension<-as.mpse(ab_hypertension)
phy<-as.phyloseq(ab_hypertension)
tmp<-prop.table(as.matrix(phy@otu_table),2)|>as.data.frame()
rownames(hypertension_meta)<-hypertension_meta$sample_id
hypertension_meta<-hypertension_meta[colnames(tmp),]
batch_rm<-adjust_batch(as.data.frame(tmp),
                       batch="study_name",
                       data=hypertension_meta)$feature_abd_adj
phy@otu_table<-otu_table(batch_rm,taxa_are_rows = T)
mm_lefse <- run_lefse(
  phy,taxa_rank="Genus",
  wilcoxon_cutoff = 0.05,
  group = "disease",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 2
)
dif_genus2<-mm_lefse@marker_table
write.csv(as.matrix(dif_genus2),file = "data/metagenome/dif_genus.csv",quote = F)
dif_genus2$feature<-gsub("g__|_f__.*|_o__.*|_p__.*","",dif_genus2$feature)
dif_genus_p2<-draw_lefse(dif_genus2)

## 16s ----
otu<-read.table("data/16s/microbio_selected.otus",row.names = 1,header = T)
meta_16s<-read.csv("data/16S/SraRunTable.csv",header = T)
rownames(meta_16s)<-meta_16s$Sample.Name
meta_16s$group<-ifelse(meta_16s$diastolic_bp>=90|meta_16s$Systolic_BP>=140,"hypertension","healthy")

meta_16s<-meta_16s[meta_16s$group!="median",]
taxa<-read.table("data/16s/microbio_selected.taxonomy",row.names = 1,header = T)
taxatmp<-strsplit(taxa$Taxonomy,";")
taxatmp<-as.data.frame(taxatmp)
taxa<-t(taxatmp)
rownames(taxa)<-colnames(otu)
colnames(taxa)<-c("Kingdom", "Phylum","Class","Order","Family" ,"Genus","Species")
taxa[,1]<-gsub("k__","",taxa[,1])
taxa[,2]<-gsub("p__","",taxa[,2])
taxa[,3]<-gsub("c__","",taxa[,3])
taxa[,4]<-gsub("o__","",taxa[,4])
taxa[,5]<-gsub("f__","",taxa[,5])
taxa[,6]<-gsub("g__","",taxa[,6])
taxa[,7]<-gsub("s__","",taxa[,7])
taxa[taxa=="unclassified"]<-NA

test_mpse<-mp_import_dada2(seqtab = otu,taxatab = taxa,sampleda = meta_16s)

phy<-as.phyloseq(test_mpse)
mm_lefse <- run_lefse(
  phy,
  wilcoxon_cutoff = 0.05,
  group = "group",
  kw_cutoff = 0.05,
  multigrp_strat = TRUE,
  lda_cutoff = 2,
  taxa_rank = "Genus"
)
dif_genus3 <- mm_lefse@marker_table
write.csv(as.matrix(dif_genus3),"data/16S/dif_genus.csv",quote = F)
dif_genus3$feature<-gsub("g__|_f__.*|_o__.*","",dif_genus3$feature)
dif_genus_p3<-draw_lefse(dif_genus3)


p1_up<-dif_genus_p1[[1]]$feature[dif_genus_p1[[1]]$lda>0]
p1_down<-dif_genus_p1[[1]]$feature[dif_genus_p1[[1]]$lda<0]
p2_up<-dif_genus_p2[[1]]$feature[dif_genus_p2[[1]]$lda>0]
p2_down<-dif_genus_p2[[1]]$feature[dif_genus_p2[[1]]$lda<0]
p3_up<-dif_genus_p3[[1]]$feature[dif_genus_p3[[1]]$lda>0]
p3_down<-dif_genus_p3[[1]]$feature[dif_genus_p3[[1]]$lda<0]
# gutMdisorder ----

gutm<-read.csv("data/gutMdisorder/Data export.csv")

gutm$Gut.Microbiota..ID.<-gsub(" \\(.*","",gutm$Gut.Microbiota..ID.)

gutm_up<-gutm$Gut.Microbiota..ID.[gutm$Alteration=="increase"]
gutm_down<-gutm$Gut.Microbiota..ID.[gutm$Alteration=="decrease"]

venup<-list(CPH16S=p1_up,
            CHMG=p2_up,
     CH16S=p3_up,
     GutMDisorder=gutm_up)

vendown<-list(CPH16S=p1_down,
            CHMG=p2_down,
            CH16S=p3_down,
            GutMDisorder=gutm_down)
pd1<-upset(fromList(venup), nsets =4, sets =names(vendown),
      keep.order =T, number.angles =0,
      line.size =1.5, mb.ratio =c(0.6,0.4),  # 主条形图和矩阵图的比例
      order.by =c("degree"), 
      decreasing = F,
      
      main.bar.color = c(rep("#000000",4) ,"#d56763","#fcd2a1") , #上方y轴柱状图颜色
      matrix.color = "#000000",             # 矩阵点颜色  
      sets.bar.color = "#4e40ae",           # 集合柱状图颜色
      point.size = 5,
      queries = list(
        list(query = intersects, params = list(names(venup[c(1,3,4)])),  color = "#d56763", active = TRUE),
        list(query = intersects, params = list(names(venup[c(1:3)])),  color = "#fcd2a1",  active = TRUE)
      ),
      text.scale =c( 2, 1.5,1.8, 2,2,2.5),
      # 标签
      sets.x.label ="",
      #mainbar.y.label ="n= 1,427 significant cell-specific\nreceptor-ligand interactions (p < 0.05)",
      
      # 边距
      set.metadata =NULL
)

pe1<-upset(fromList(vendown), nsets =4, sets =names(vendown),
      keep.order =T, number.angles =0,
      line.size =1.5, mb.ratio =c(0.6,0.4),  # 主条形图和矩阵图的比例
      order.by =c("degree"), 
      decreasing = F,
      
      main.bar.color = c(rep("#000000",4) ,"#d56763","#fcd2a1","#477b80","#2aa080","#a760ae") , #上方y轴柱状图颜色
      matrix.color = "#000000",             # 矩阵点颜色  
      sets.bar.color = "#4e40ae",           # 集合柱状图颜色
      point.size = 5,
      queries = list(
        list(query = intersects, params = list(names(vendown[c(2,4)])),  color = "#d56763", active = TRUE),
        list(query = intersects, params = list(names(vendown[c(2,3)])),  color = "#fcd2a1",  active = TRUE),
        list(query = intersects, params = list(names(vendown[c(1,3)])),  color = "#477b80",  active = TRUE),
        list(query = intersects, params = list(names(vendown[c(1,2)])),  color = "#2aa080",  active = TRUE),
        list(query = intersects, params = list(names(vendown[c(1:3)])),  color = "#a760ae",  active = TRUE) 
      ),
      
      # text.scale参数说明：
      # • 第1个值：Y轴标题
      # • 第2个值：Y轴刻度标签
      # • 第3个值：X轴标题
      # • 第4个值：X轴刻度标签
      # • 第5个值：集合名字大小
      # • 第6个值：柱体s上的数字
      text.scale =c( 2, 1.5,1.8, 2,2,2.5),
      # 标签
      sets.x.label ="",
      #mainbar.y.label ="n= 1,427 significant cell-specific\nreceptor-ligand interactions (p < 0.05)",
      
      # 边距
      set.metadata =NULL
      )

pa<-dif_genus_p1[[2]]+ggtitle("A")+
  theme(legend.position = "none")+labs(x="")
pb<-dif_genus_p2[[2]]+ggtitle("B")+
  theme(legend.position = "bottom")
pc<-dif_genus_p3[[2]]+ggtitle("C")+
  theme(legend.position = "none")

tmp1<-grid.arrange(grid.arrange(pa,pc,nrow=2,heights=c(0.55,0.45)),pb,nrow=1)

pd<-grid.grabExpr(print(pd1))
pe<-grid.grabExpr(print(pe1))
tmp2<-cowplot::plot_grid(pd, pe, ncol = 1, labels = c("D", "E"), label_size = 18)
dev.off()
pdf("result/dif_species/dif_genus_all(Fig1).pdf",width = 14,height = 10)
grid.arrange(tmp1,tmp2,nrow=1,widths=c(6,3.5))
dev.off()

combs<-combn(names(vendown),2)
inter_up<-apply(combs,2,function(x){
  Reduce(intersect,venup[x])
})|>unlist()|>unique()
inter_down<-apply(combs,2,function(x){
  Reduce(intersect,vendown[x])
})|>unlist()|>unique()
dif_res<-data.frame(genus=c(inter_up,inter_down),
           enrich_group=rep(c("hypertension","healthy"),
                            c(length(inter_up),length(inter_down))))
write.csv(dif_res,"result/dif_species/dif_genus_intersect.csv",row.names = F,quote = F)
