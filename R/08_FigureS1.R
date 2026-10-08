# 07_FigureS1.R
source("R/00_setup.R")
d <- read_csv(file.path(D_PANEL,"Supplementary","FigureS1",
                        "FigS1_day0_extended_alpha_data.csv"),show_col_types=FALSE)
metrics <- c("Bacterial_richness","Bacterial_Chao1","Bacterial_PD",
             "Fungal_richness","Fungal_Chao1","Fungal_PD",
             "KO_Shannon","COG_Shannon")
z <- d %>% select(SampleID,Treatment,any_of(metrics)) %>%
  pivot_longer(-c(SampleID,Treatment),names_to="Metric",values_to="Value") %>%
  drop_na()
p <- ggplot(z,aes(Treatment,Value,colour=Treatment))+
  geom_jitter(width=.08,size=1.3)+
  stat_summary(fun=mean,geom="crossbar",width=.45,linewidth=.35,colour="black")+
  facet_wrap(~Metric,scales="free_y",ncol=4)+
  scale_colour_manual(values=treatment_cols)+
  labs(x=NULL,y=NULL)+theme_pub()+theme(legend.position="none")
save_pdf(p,"FigureS1.pdf",width=7.2,height=4.2)
