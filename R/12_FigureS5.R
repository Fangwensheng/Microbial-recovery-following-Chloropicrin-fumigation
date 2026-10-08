# 11_FigureS5.R
source("R/00_setup.R")
driver_file <- file.path(D_REC,"Fig5_driver_master_reclassified_Rrerun.csv")
if(!file.exists(driver_file)) driver_file <- file.path(D_REC,"Fig5_driver_master_reclassified_v2.csv")
driver <- read_csv(driver_file,show_col_types=FALSE)
day0 <- read_csv(file.path(D_RAW,"DESeq2_day0_all_KO.csv"),show_col_types=FALSE)
rec <- read_csv(file.path(D_RAW,"DESeq2_KO_recovery_descriptors.csv"),show_col_types=FALSE)
resolved <- driver$KO

# a: MAG-resolved selection bias at day 0
a0 <- day0 %>% filter(!is.na(padj),padj<0.10,abs(log2FoldChange)>=0.5) %>%
  transmute(KO,Initial_abs_log2FC=abs(log2FoldChange),
            MAG_resolved=if_else(KO %in% resolved,"Resolved","Unresolved"))
a <- ggplot(a0,aes(MAG_resolved,Initial_abs_log2FC,fill=MAG_resolved))+
  geom_boxplot(outlier.shape=NA,width=.6)+
  scale_fill_manual(values=c(Resolved="#4C78A8",Unresolved="grey70"))+
  theme_pub()+theme(legend.position="none")+
  labs(x=NULL,y=expression("|day-0 log"[2]*"FC|"))

# b: recovery residuals in resolved vs unresolved KOs
b0 <- rec %>% mutate(MAG_resolved=if_else(KO %in% resolved,"Resolved","Unresolved"),
                     Context=recode(Mode,Natural="Unplanted",Plant="Tomato-planted"))
b <- ggplot(b0,aes(MAG_resolved,Time_integrated_absolute_residual_7to70d,fill=MAG_resolved))+
  geom_boxplot(outlier.shape=NA,width=.6)+facet_wrap(~Context)+
  scale_fill_manual(values=c(Resolved="#4C78A8",Unresolved="grey70"))+
  theme_pub()+theme(legend.position="none")+
  labs(x=NULL,y="Cumulative residual")

# c: correlations among redundancy metrics / major static predictors
cols <- c("Provider_richness","Effective_provider_diversity",
          "Provider_phylogenetic_breadth","Weighted_KO_copy",
          "Weighted_GC_percent","Genome_size_Mb","Weighted_rrn_copy",
          "Weighted_CAZyme_density")
cm <- cor(driver[,cols],use="pairwise.complete.obs",method="spearman")
cm_long <- as.data.frame(as.table(cm)) %>% rename(X=Var1,Y=Var2,rho=Freq)
c <- ggplot(cm_long,aes(X,Y,fill=rho))+geom_tile()+
  scale_fill_gradient2(low="#3B6FB6",mid="white",high="#C94C4C",midpoint=0,limits=c(-1,1))+
  theme_pub()+theme(axis.text.x=element_text(angle=45,hjust=1))+
  labs(x=NULL,y=NULL,fill=expression(rho))

# d: adjusted redundancy associations
pp_file <- file.path(D_REC,"Fig5_partial_spearman_Rrerun.csv")
if(!file.exists(pp_file)) pp_file <- file.path(D_REC,"Fig5_partial_spearman_v2.csv")
pp <- read_csv(pp_file,show_col_types=FALSE)
d <- ggplot(pp,aes(partial_rho,Predictor,colour=Context,shape=Context))+
  geom_vline(xintercept=0,linetype=2,colour="grey60")+geom_point(size=2)+
  scale_colour_manual(values=c(Unplanted="#4C6F8A",`Tomato-planted`="#2A9D8F"))+
  theme_pub()+labs(x="Adjusted partial Spearman \u03c1",y=NULL)

save_pdf((a|b)/(c|d)+plot_annotation(tag_levels="a"),"FigureS5_R.pdf",width=7.2,height=6.8)
