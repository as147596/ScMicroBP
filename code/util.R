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
    scale_fill_manual(values=c("#F9D5B2","#9FC3E2"))
  list(lefse,p)
}
