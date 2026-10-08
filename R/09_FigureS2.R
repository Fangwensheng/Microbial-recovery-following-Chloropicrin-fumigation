# 08_FigureS2_v2.R

source("R/00_setup_Fig2_S2_v2.R")

a <- read_csv(file.path(D,"FigS2_alpha_diversity_raw.csv"),show_col_types=FALSE) %>%
  mutate(Treatment=factor(Treatment,levels=c("CK","CP","TCK","TCP")))
pc <- read_csv(file.path(D,"FigS2_PCoA_coordinates.csv"),show_col_types=FALSE) %>%
  filter(Day>=7) %>%
  mutate(Treatment=factor(Treatment,levels=c("CK","CP","TCK","TCP")))

metric_info <- tribble(
 ~Metric,~Label,~Ylab,
 "Bac_Richness","Bacterial richness","Observed ASVs",
 "Fun_Richness","Fungal richness","Observed ASVs",
 "Bac_Chao1","Bacterial Chao1","Chao1",
 "Fun_Chao1","Fungal Chao1","Chao1",
 "Bac_PD","Bacterial Faith's PD","Faith's PD",
 "Fun_PD","Fungal Faith's PD","Faith's PD"
)

al <- a %>%
  pivot_longer(cols=all_of(metric_info$Metric),names_to="Metric",values_to="Value") %>%
  group_by(Treatment,Day,Metric) %>%
  summarise(mean=mean(Value),se=sd(Value)/sqrt(n()),.groups="drop") %>%
  left_join(metric_info,by="Metric")

make_alpha <- function(metric){
  z <- al %>% filter(Metric==metric)
  info <- metric_info %>% filter(Metric==metric)
  ggplot(z,aes(Day,mean,colour=Treatment,shape=Treatment,
               linetype=Treatment,group=Treatment))+
    geom_line(linewidth=.58)+geom_point(size=1.55)+
    geom_errorbar(aes(ymin=mean-se,ymax=mean+se),width=1.6,linewidth=.32)+
    scale_colour_manual(values=treatment_cols,drop=FALSE)+
    scale_shape_manual(values=treatment_shapes,drop=FALSE)+
    scale_linetype_manual(values=treatment_ltypes,drop=FALSE)+
    scale_x_continuous(breaks=c(7,14,28,42,56,70))+
    labs(x="Days after fumigation",y=info$Ylab)+theme_nc()+
    annotate("text",x=7,y=Inf,label=info$Label,hjust=-0.02,vjust=1.35,size=2.7)
}

a1<-make_alpha("Bac_Richness"); a2<-make_alpha("Fun_Richness")
a3<-make_alpha("Bac_Chao1");    a4<-make_alpha("Fun_Chao1")
a5<-make_alpha("Bac_PD");       a6<-make_alpha("Fun_PD")

cent <- pc %>%
  group_by(Profile,Treatment,Day) %>%
  summarise(PCoA1=mean(PCoA1),PCoA2=mean(PCoA2),.groups="drop")
axisinfo <- pc %>%
  group_by(Profile) %>%
  summarise(xpct=first(PCoA1_percent),ypct=first(PCoA2_percent),.groups="drop")

make_pcoa <- function(profile){
  z<-cent %>% filter(Profile==profile)
  ai<-axisinfo %>% filter(Profile==profile)
  ggplot(z,aes(PCoA1,PCoA2,colour=Treatment,shape=Treatment,
               linetype=Treatment,group=Treatment))+
    geom_path(linewidth=.58,
              arrow=grid::arrow(length=grid::unit(0.055,"inches"),type="closed"))+
    geom_point(size=1.55)+
    scale_colour_manual(values=treatment_cols,drop=FALSE)+
    scale_shape_manual(values=treatment_shapes,drop=FALSE)+
    scale_linetype_manual(values=treatment_ltypes,drop=FALSE)+
    labs(x=sprintf("PCoA1 (%.1f%%)",ai$xpct),
         y=sprintf("PCoA2 (%.1f%%)",ai$ypct))+
    theme_nc()+
    annotate("text",x=-Inf,y=Inf,label=paste0(profile," Aitchison PCoA"),
             hjust=-0.05,vjust=1.35,size=2.7)+
    annotate("text",x=-Inf,y=Inf,label="Descriptive centroid trajectory",
             hjust=-0.05,vjust=3.1,size=2.05,colour="grey35")
}

p7<-make_pcoa("Bacteria"); p8<-make_pcoa("Fungi")
p <- ((a1|a2)/(a3|a4)/(a5|a6)/(p7|p8))+
  plot_annotation(tag_levels="a",
                  theme=theme(plot.tag=element_text(face="bold",size=9)))+
  plot_layout(guides="collect") &
  theme(legend.position="top")

save_both(p,"FigureS2",width=7.2,height=10.0)
