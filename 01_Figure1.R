# 01_Figure1.R
source("R/00_setup.R")

raw <- read_csv(file.path(D_PANEL,"Main_Figures","Figure1","Fig1a_raw_scalar_data.csv"),
                show_col_types=FALSE)
eff <- read_csv(file.path(D_PANEL,"Main_Figures","Figure1","Fig1b_scalar_effects_and_statistics.csv"),
                show_col_types=FALSE)
pc <- read_csv(file.path(D_PANEL,"Main_Figures","Figure1","Fig1c_PCoA_coordinates.csv"),
               show_col_types=FALSE)

# a: day-0 scalar observations. qPCR abundance is log10-transformed only for plotting.
a <- raw %>%
  mutate(Plot_value=if_else(str_detect(Dimension,"gene abundance"),
                            log10(Value),Value),
         Dimension=factor(Dimension,levels=c(
           "16S rRNA gene abundance","ITS gene abundance",
           "Bacterial Shannon","Fungal Shannon","KO richness","COG richness"))) %>%
  ggplot(aes(Treatment,Plot_value,colour=Treatment))+
  geom_jitter(width=.08,height=0,size=1.5)+
  stat_summary(fun=mean,geom="crossbar",width=.45,linewidth=.35,colour="black")+
  facet_wrap(~Dimension,scales="free_y",nrow=2)+
  scale_colour_manual(values=treatment_cols)+
  labs(x=NULL,y=NULL)+theme_pub()+theme(legend.position="none")

# b: effect-size overview (log2 CP/CK for all scalar metrics).
b <- eff %>%
  mutate(Profile_or_metric=factor(Profile_or_metric,
                                 levels=rev(Profile_or_metric))) %>%
  ggplot(aes(Effect_or_F,Profile_or_metric))+
  geom_vline(xintercept=0,linetype=2,colour="grey60")+
  geom_point(aes(fill=FDR_or_PERMDISP_P<0.05),shape=21,size=2.1)+
  scale_fill_manual(values=c(`TRUE`=COL_FUM,`FALSE`="white"))+
  labs(x=expression(log[2]*"(CP / CK)"),y=NULL)+theme_pub()+
  theme(legend.position="none")

# c: PCoA coordinates. PCoA axis sign is arbitrary; flipping an axis has no
# biological consequence.
c <- pc %>%
  ggplot(aes(PCoA1,PCoA2,colour=Treatment,shape=Treatment))+
  geom_point(size=2)+
  facet_wrap(~Profile,scales="free",nrow=2,
             labeller=labeller(Profile=function(x)x))+
  scale_colour_manual(values=treatment_cols)+
  scale_shape_manual(values=treatment_shapes)+
  labs(x="PCoA1",y="PCoA2")+theme_pub()

p <- a / (b | c) + plot_annotation(tag_levels="a")
save_pdf(p,"Figure1.pdf",width=7.2,height=7.4)
