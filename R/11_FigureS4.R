# 10_FigureS4.R
source("R/00_setup.R")
d <- read_csv(file.path(D_PANEL,"Supplementary","FigureS4",
                        "FigS4_DESeq2_and_clustering_sensitivity.csv"),show_col_types=FALSE)
rec <- read_csv(file.path(D_RAW,"DESeq2_KO_recovery_descriptors.csv"),show_col_types=FALSE)
primary <- read_csv(file.path(D_RAW,"DESeq2_disturbed_KO_primary.csv"),show_col_types=FALSE)
hc <- read_csv(file.path(D_RAW,"DESeq2_day0_all_KO.csv"),show_col_types=FALSE) %>%
  filter(!is.na(padj),padj<0.05,abs(log2FoldChange)>=1)

v <- d %>% filter(Section=="DESeq2 day0 sensitivity") %>%
  mutate(mlog10p=-log10(pmax(pvalue,1e-300)),
         sig=Primary_disturbed)
p1 <- ggplot(v,aes(log2FoldChange,mlog10p,colour=sig))+
  geom_point(alpha=.4,size=.7)+geom_vline(xintercept=c(-.5,.5),linetype=2)+
  scale_colour_manual(values=c(`TRUE`=COL_FUM,`FALSE`="grey75"))+
  theme_pub()+theme(legend.position="none")+
  labs(x=expression("Day-0 log"[2]*"(CP/CK)"),y=expression(-log[10]*P))

sil <- d %>% filter(Section=="Silhouette selection")
p2 <- ggplot(sil,aes(k,Silhouette))+geom_line()+geom_point()+
  theme_pub()+labs(x="k",y="Silhouette coefficient")

sensitivity <- bind_rows(
  rec %>% filter(KO %in% primary$KO) %>% mutate(Set="Primary"),
  rec %>% filter(KO %in% hc$KO) %>% mutate(Set="High-confidence")
) %>%
  group_by(Set,Mode) %>% summarise(Median=median(Time_integrated_absolute_residual_7to70d),
                                   N=n(),.groups="drop")
p3 <- ggplot(sensitivity,aes(Set,Median,fill=Mode))+geom_col(position="dodge")+
  theme_pub()+labs(x=NULL,y="Median cumulative residual")

end <- rec %>% select(KO,Mode,Absolute_residual_at_70d) %>%
  pivot_wider(names_from=Mode,values_from=Absolute_residual_at_70d)
p4 <- ggplot(end,aes(Natural,Plant))+geom_abline(slope=1,intercept=0,linetype=2)+
  geom_point(alpha=.2,size=.7)+theme_pub()+
  labs(x="Unplanted 70-d residual",y="Tomato-planted 70-d residual")
save_pdf((p1|p2)/(p3|p4)+plot_annotation(tag_levels="a"),"FigureS4_R.pdf",width=7.2,height=6.8)
