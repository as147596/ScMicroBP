plot_graph<-function(net,nodes){
  graph<-graph_from_data_frame(net,vertices = nodes)
  V(graph)$degree<-igraph::degree(graph)
  nodes$degree<-igraph::degree(graph)
  nodes<-nodes[order(nodes$cluster),]
  graph<-graph_from_data_frame(net,vertices = nodes)
  cluster<-props$clustering$clust1
  V(graph)$enrich_group<-factor(V(graph)$enrich_group,levels = c("healthy","hypertension","Notsignificant"))
  V(graph)$lab<-ifelse(V(graph)$enrich_group!="Notsignificant",names(V(graph)),NA)
  V(graph)$modularity3<-V(graph)$cluster
  V(graph)$Modularity<-V(graph)$cluster
  V(graph)$Degree<-igraph::degree(graph)
  graph_tb<-tidygraph::as_tbl_graph(graph)
  p1<-ggNetView(
    graph_obj = graph_tb,
    layout = "gephi",
    layout.module = "random",
    pointsize = c(1, 5),
    group.by = "Modularity",
    center = F,
    jitter = F,
    mapping_line = F,
    shrink = 0.9,
    linealpha = 0.2,
    linecolor = "#d9d9d9"
  ) 
  graph<-graph_from_data_frame(net,vertices = p1[["layers"]][[2]][["data"]])
  
  ggraph(graph,
         layout = "manual",
         x = x,
         y = y) +
    geom_edge_link0(aes(col = association),width=0.2) +
    geom_node_point(aes(shape = enrich_group,fill = cluster,color=cluster,  size = degree)) +
    scale_shape_manual(values = c(Notsignificant=1,healthy=19,hypertension=17))+
    geom_mark_hull( 
      aes(x, y, group = cluster, fill = cluster),
      concavity = 4,
      expand = unit(2, "mm"),
      alpha = 0.1
    ) +
    geom_node_text(aes(label = lab),repel = T)+
    scale_color_manual(values = c("orange3","orchid4","red","orangered4",
                                  "palegreen4","royalblue4","orange","linen",
                                  "skyblue4","thistle","turquoise4","lightsteelblue4","mediumseagreen")
    ) +
    scale_fill_manual(values = c("orange3","orchid4","red","orangered4",
                                 "palegreen4","royalblue4","orange","linen",
                                 "skyblue4","thistle","turquoise4","lightsteelblue4","mediumseagreen")
    ) +
    scale_edge_colour_manual(
      values = c("aquamarine3","mediumpurple"))+
    guides(fill = guide_legend(ncol = 1))+
    theme_graph(base_family = "sans")+
    #theme(legend.box = 'horizontal',
    #      legend.box.just = 'top')+
    theme(plot.title = element_text(hjust = -0.02,vjust = 4.5),
          title = element_text(size = 16,face="bold"))
  
}

draw_lefse<-function(lefse){
  lefse<-lefse[lefse$feature!="un",]
  lefse<-lefse[!duplicated(lefse$feature),]
  lefse$lda<-ifelse(lefse$enrich_group=="healthy",
                    -lefse$ef_lda,
                    lefse$ef_lda)
  marker_table<-lefse
  p<-ggplot(marker_table,aes(x=lda,y=reorder(feature,lda),fill=enrich_group))+
    geom_col()+
    theme_test()+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_text(size = 16,face = "bold"),
          legend.text = element_text(size = 14),
          axis.title = element_text(size=14,face = "bold"),
          axis.text = element_text(size=10),
          plot.margin = margin(t=5,b=5,l=10,r=10),
          legend.position = "right",
          plot.title = element_text(hjust = -0.065,vjust = 0.5),
          title = element_text(size = 16,face="bold"))+
    geom_text(data = marker_table[marker_table$enrich_group == "healthy",], aes(y = feature, x = 0.1, label = feature),
              hjust = 0, size = 3.5,fontface = "italic") +
    geom_text(data = marker_table[marker_table$enrich_group == "hypertension",], aes(y = feature, x = -0.1, label = feature),
              hjust = 1, size = 3.5,fontface="italic")+expand_limits(x=c(-6,5))+
    labs(x="LDA score",y="")+
    scale_fill_manual(values=c("#9FC3E2","#F9D5B2"))
  list(lefse,p)
}
cell_DEG_enrich<-function(srt,group.by="group",top=20){
  DEG<-FindMarkers(srt,ident.1 = "Hypertension",ident.2 = "Control",
                   group.by = group.by,logfc.threshold=0,min.pct=0)
  DEG$Change=ifelse(DEG$p_val_adj<0.05&abs(DEG$avg_log2FC)>0.58,ifelse(DEG$avg_log2FC>0.58,"Up","Down"),"Not")
  colors<- c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25")
  p<-ggplot(DEG,aes(x=avg_log2FC,y= -log10(p_val_adj)))+
    geom_point(aes(color = Change, size = -log10(p_val_adj)), alpha =0.5) +
    scale_color_manual(values = c(Down="#39489f",Not="grey80",Up="#b81f25")) +
    ggthemes::theme_clean()+
    geom_hline(yintercept = -log10(0.05),linetype="dashed",size=1)+
    geom_vline(xintercept = c(-0.58,0.58),linetype=2,size=1)+
    theme(legend.key = element_rect(fill = 'transparent'),
          legend.background = element_rect(fill = 'transparent'),
          axis.title = element_text(size = 16),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.margin = ggplot2::margin(r=10,b=10,t=10,l=10))+
    theme(axis.text = element_text(color = "black",size = 12),
          axis.title = element_text(size = 14,face = "bold"),
          plot.title.position = "plot",
          plot.title = element_text(size=14,vjust = -2,hjust = -0.01,face = "bold"),
          strip.text = element_text(face = "bold",size = 16),
          legend.key.size = unit(0.5,'cm'), 
          legend.title = element_text(size = 14,face = "bold"),
          plot.background = element_blank(),
          legend.text = element_text(size = 12),
          legend.background = element_blank())
  genes<-bitr(rownames(DEG)[DEG$Change!="Not"],fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
  kk<-enrichKEGG(genes$ENTREZID)
  kk_df<-kk@result[kk@result$pvalue<0.05,]|>as.data.frame()
  if(nrow(kk_df)>top){
    kk_df<-kk_df[1:top,]
  }
  kk_df$Ratio=as.numeric(gsub("/.*","",kk_df$GeneRatio))/as.numeric(gsub(".*/","",kk_df$GeneRatio))
  ep<-ggplot(kk_df,aes(Ratio,reorder(Description,Ratio)))+
    geom_point(aes(color=pvalue,size=Count))+
    scale_colour_distiller(palette="OrRd")+
    theme_classic()+
    labs(y="")+
    theme(legend.key = element_rect(fill = 'transparent'),
          legend.background = element_rect(fill = 'transparent'),
          axis.title = element_text(size = 16),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.margin = ggplot2::margin(r=10,b=10,t=10,l=10))+
    theme(axis.text = element_text(color = "black",size = 12),
          axis.title = element_text(size = 14,face = "bold"),
          plot.title.position = "plot",
          plot.title = element_text(size=14,vjust = -2,hjust = -0.01,face = "bold"),
          strip.text = element_text(face = "bold",size = 16),
          legend.key.size = unit(0.5,'cm'), 
          legend.title = element_text(size = 14,face = "bold"),
          plot.background = element_blank(),
          legend.text = element_text(size = 12),
          legend.background = element_blank())
  return(list(p,ep))
}
DEP<-function(str,celltype,aucMatrix){
  bc<-str$subcelltype==celltype
  list <- str$group[bc]|>factor(levels = c("Control", "Hypertension"))
  list <- model.matrix(~factor(list)+0)  #把group设置成一个model matrix
  colnames(list) <- c("Control", "Hypertension")
  df.fit <- lmFit(t(aucMatrix[colnames(str)[bc],]), list)
  df.matrix <- makeContrasts(Hypertension - Control , levels = list)
  fit <- contrasts.fit(df.fit, df.matrix)
  fit <- eBayes(fit)
  DEP <- topTable(fit,n = Inf, adjust = "fdr")
  DEP$lab<-ifelse(rownames(DEP)%in%c(
    "HALLMARK_INFLAMMATORY_RESPONSE","HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_IL2_STAT5_SIGNALING"
  ),rownames(DEP),"")
  DEP$lab[DEP$adj.P.Val>0.05]<-""
  
  colors<- c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25")
  p<-ggplot(DEP,aes(x=logFC,y= -log10(adj.P.Val)))+
    geom_point(aes(color = logFC, size = -log10(adj.P.Val)), alpha =0.5) +
    scale_color_gradientn(colors = colors) +
    ggthemes::theme_clean()+
    geom_text_repel(aes(label = lab),max.overlaps = Inf,box.padding = 1)+
    geom_hline(yintercept = -log10(0.05),linetype="dashed",size=1)+
    geom_vline(xintercept = 0,linetype=2,size=1)+
    theme(legend.key = element_rect(fill = 'transparent'),
          legend.background = element_rect(fill = 'transparent'),
          axis.title = element_text(size = 16),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.margin = ggplot2::margin(r=10,b=10,t=10,l=10))+
    theme(axis.text = element_text(color = "black",size = 12),
          axis.title = element_text(size = 14,face = "bold"),
          plot.title.position = "plot",
          plot.title = element_text(size=14,vjust = -2,hjust = -0.01,face = "bold"),
          strip.text = element_text(face = "bold",size = 16),
          legend.key.size = unit(0.5,'cm'), 
          legend.title = element_text(size = 14,face = "bold"),
          plot.background = element_blank(),
          legend.text = element_text(size = 12),
          legend.background = element_blank())
  list(DEP,p)
}

pathway_cor<-function(srt,sp="Bacteroides",cor_res,top=3,sel=NULL,kegg_aucMatrix=kegg_aucMatrix,
                      ylab="Mean score of KEGG pathways"){
  if(is.null(sel)){ind<-order(abs(cor_res),decreasing = T)[1:top]}
  else{
    ind<-sel
  }
  kegg_aucMatrix1<-kegg_aucMatrix[,ind,drop=F]
  tmp<-srt@meta.data
  kegg_aucMatrix1[,sp]<-case_when(tmp[,sp]<0~0,
                                  tmp[,sp]>0&tmp[,sp]<3~2,
                                  tmp[,sp]>=3&tmp[,sp]<5~4,
                                  tmp[,sp]>=5&tmp[,sp]<7~6,
                                  tmp[,sp]>=7&tmp[,sp]<9~8,
                                  tmp[,sp]>=9~10,
  )
  kegg_aucMatrix1<-aggregate(as.formula(paste0(".~",sp)),kegg_aucMatrix1,mean)
  rownames(kegg_aucMatrix1)<-kegg_aucMatrix1[,1]
  kegg_aucMatrix1<-kegg_aucMatrix1[,-1]
  kegg_aucMatrix1_long<-reshape2::melt(as.matrix(kegg_aucMatrix1))
  kegg_aucMatrix1_long$genus<-gsub("_BPS","",sp)
  p<-ggplot(kegg_aucMatrix1_long,aes(x=Var1))+
    geom_line(aes(y=value,color=Var2),linewidth = 1.5)+
    geom_point(aes(y=value,color=Var2),size=4)+
    labs(x=paste0(gsub("_BPS","",sp)," ","BPS (quintile bins)"),
         y=ylab)+
    theme_classic()+
    scale_color_brewer(name = "",palette = "Set3")+
    scale_x_continuous(breaks = seq(0,10,2))
  list(p,kegg_aucMatrix1_long)
}

multi_vol<-function(x1,x2,x3,x4,cells=c("CD4 Naive","CD4 TCM","CD8 Naive","CD8 TCM")){
  addy=150
  addx=0.04
  x1$y=-log10(x1$adj.P.Val)+addy
  x1$x=x1$logFC
  x2$y=-log10(x2$adj.P.Val)
  x2$x=x2$logFC
  x1$side<-x2$side<-"left"
  x3$y=-log10(x3$adj.P.Val)+addy
  x3$x=-x3$logFC+addx
  x4$y=-log10(x4$adj.P.Val)
  x4$x=-x4$logFC+addx
  x3$side<-x4$side<-"right"
  
  lab<-rbind(x1[x1$lab!="",],
             x2[x2$lab!="",],
             x3[x3$lab!="",],
             x4[x4$lab!="",])
  y=seq(50,350,length=5)
  xe<-rep(0.012,5)
  pos<-data.frame(pathway=unique(lab$lab),
                  xend=xe,
                  yp=y)
  lab$xend=pos$xend[match(lab$lab,pos$pathway)]
  lab$yp=pos$yp[match(lab$lab,pos$pathway)]
  tit<-data.frame(x=-0.001,0.001)
  colors<- c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25")
  ggplot()+
    geom_hline(yintercept = 150)+
    geom_point(data = x1,mapping = aes(logFC, -log10(adj.P.Val)+addy,
                                       size = -log10(adj.P.Val),colour = logFC))+
    geom_point(data = x2,mapping = aes(logFC, -log10(adj.P.Val),
                                       size = -log10(adj.P.Val),colour = logFC))+
    geom_point(data = x3,mapping = aes(-logFC+addx, -log10(adj.P.Val)+addy,
                                       size = -log10(adj.P.Val),colour = logFC))+
    geom_point(data = x4,mapping = aes(-logFC+addx, -log10(adj.P.Val),
                                       size = -log10(adj.P.Val),colour = logFC))+
    annotate("rect",xmin = xe[1], xmax = xe[1]+0.016, ymin = -Inf, ymax = Inf, 
             fill = "#8dd3c7", alpha = 0.3)+
    geom_text(data = lab,mapping=aes(0.02,yp,label = lab),size=3)+
    geom_point(data = lab,mapping=aes(xend,yp),size=5,color="#a38466")+
    geom_point(data = lab,mapping=aes(xend+0.016,yp),size=5,color="#a38466")+
    geom_segment(data = lab[lab$side=="left",],mapping =aes(x=x,y=y,xend = xend,yend=yp),linetype = 2,color="grey50")+
    geom_segment(data = lab[lab$side=="right",],mapping =aes(x=x,y=y,xend = xend+0.016,yend=yp),linetype = 2,color="grey50")+
    scale_color_gradientn(colors = colors) +
    ggthemes::theme_clean()+
    geom_hline(yintercept = c(-log10(0.05),-log10(0.05)+150),linetype=2)+
    scale_x_continuous(breaks = c(seq(-0.0025,0.0115,0.005),
                                  c(0.029,0.035,0.04,0.044)),
                       labels = c(seq(-0.0025,0.0115,0.005),
                                  c(0.011,0.005,0,-0.009)))+
    scale_y_continuous(breaks = c(seq(0,140,40),seq(150,400,40)),
                       labels = c(seq(0,140,40),seq(150,400,40)-150))+
    theme(legend.key = element_rect(fill = 'transparent'),
          legend.background = element_rect(fill = 'transparent'),
          axis.title = element_text(size = 16),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.margin = ggplot2::margin(r=10,b=10,t=10,l=10))+
    theme(axis.text = element_text(color = "black",size = 12),
          axis.title = element_text(size = 14,face = "bold"),
          plot.title.position = "plot",
          plot.title = element_text(size=14,vjust = -2,hjust = -0.01,face = "bold"),
          strip.text = element_text(face = "bold",size = 16),
          legend.key.size = unit(0.5,'cm'), 
          legend.title = element_text(size = 14,face = "bold"),
          plot.background = element_blank(),
          legend.text = element_text(size = 12),
          legend.background = element_blank())+
    labs(x="logFC",y="-log10(adj.P.Pvalue)")+
    annotate("text",x = c(0.0025,0.0025,0.038,0.038), y=c(400,130,400,130), label=cells)
}



double_vis<-function(res,celltype="dnT"){
  res$fdr<-p.adjust(res$pvalue,method = "BH")
  res$label<-paste0(round(res$ATE,3),"(",round(res$lower,3),"~",round(res$upper,3),")")
  res$id<-1:nrow(res)
  res$dir<-ifelse(res$ATE>0,"Positive","Negtive")
  res$FDR<-ifelse(res$fdr<0.001,"<0.001",res$fdr)
  res$outcome<-factor(res$outcome,levels = res$outcome[order(res$ATE)])
  if(celltype=="Monocyte"){
    ggplot(res,aes(y=reorder(outcome,ATE)))+
      geom_rect(data=res %>% filter(id %% 2 == 0),inherit.aes = F,
                aes(ymin=id-0.5,ymax=id+0.5),xmin=-1,xmax=Inf,
                fill="honeydew3",alpha=0.5)+
      geom_point(aes(x=ATE,y=outcome,color=dir),shape=15,size=3) +
      scale_color_manual(name="",values = c("#9FC3E2","#F9D5B2"))+
      geom_errorbarh(aes(xmin=lower,xmax=upper),height=0.2) +
      labs(x=NULL,y=NULL)+ 
      coord_cartesian(clip = "off")+
      scale_x_continuous(limits = c(-1,-0.3),
                         breaks = seq(-0.4,-0.3,0.1),
                         labels = seq(-0.4,-0.3,0.1),
                         guide=guide_axis(cap = TRUE)) +
      geom_text(data=res,
                aes(x =-0.6,y=outcome,label=label),size=4,fontface = "italic")+
      geom_text(data=res ,
                fontface = "italic",hjust=0,
                aes(x =-0.5,y=outcome,label=FDR),
                size=4,color="black")+
      geom_text(aes(x =-1,y=outcome,label=outcome),
                size=4,fontface = "italic",hjust=0) +
      theme(axis.text.y = element_blank())+
      geom_text(aes(x =-0.6,y=max(id)+1),size=4,label="ATE (95% CI)",fontface = "bold") +
      geom_text(aes(x =-1,y=max(id)+1),size=4,hjust=0,
                label="Outcome",fontface="bold") +
      geom_text(aes(x =-0.5,y=max(id)+1),size=4,
                label="FDR",fontface="bold")+
      annotate(geom="segment",x=0,xend=0,
               y=0.5,yend=4.5,color="black",linetype=2)+
      theme(axis.text.y=element_blank(),
            axis.text.x=element_text(color="black",face="bold",size=11),
            axis.ticks.y=element_blank(),
            panel.background = element_blank(),
            axis.line.x=element_line(color="black"),
            plot.margin = margin(r=0.5,l=0,b=.5,t=.5,unit="cm"),
            legend.position = "none")
  }
  if(celltype=="dnT"){
    ggplot(res,aes(y=reorder(outcome,ATE)))+
      geom_rect(data=res %>% filter(id %% 2 == 0),inherit.aes = F,
                aes(ymin=id-0.5,ymax=id+0.5),xmin=-3.5,xmax=Inf,
                fill="honeydew3",alpha=0.5)+
      geom_point(aes(x=ATE,y=outcome,color=dir),shape=15,size=3) +
      scale_color_manual(values = c("#9FC3E2","#F9D5B2"))+
      geom_errorbarh(aes(xmin=lower,xmax=upper),height=0.2) +
      labs(x=NULL,y=NULL)+ 
      coord_cartesian(clip = "off")+
      scale_x_continuous(limits = c(-3.5,0.4),
                         breaks = seq(-0.4,0.4,0.2),
                         labels = seq(-0.4,0.4,0.2),
                         guide=guide_axis(cap = TRUE)) +
      geom_text(data=res,
                aes(x =-1.35,y=outcome,label=label),size=4,fontface = "italic")+
      geom_text(data=res ,
                fontface = "italic",hjust=0,
                aes(x =-0.8,y=outcome,label=FDR),
                size=4,color="black")+
      geom_text(aes(x =-3.5,y=outcome,label=outcome),
                size=4,fontface = "italic",hjust=0) +
      theme(axis.text.y = element_blank())+
      geom_text(aes(x =-1.35,y=max(id)+1),size=4,label="ATE (95% CI)",fontface = "bold") +
      geom_text(aes(x =-3.5,y=max(id)+1),size=4,hjust=0,
                label="Outcome",fontface="bold") +
      geom_text(aes(x =-0.7,y=max(id)+1),size=4,
                label="FDR",fontface="bold")+
      annotate(geom="segment",x=0,xend=0,
               y=0.5,yend=4.5,color="black",linetype=2)+
      theme(axis.text.y=element_blank(),
            axis.text.x=element_text(color="black",face="bold",size=11),
            axis.ticks.y=element_blank(),
            panel.background = element_blank(),
            axis.line.x=element_line(color="black"),
            plot.margin = margin(r=0.5,l=0,b=.5,t=.5,unit="cm"),
            legend.position = "none")
  }
  
}
