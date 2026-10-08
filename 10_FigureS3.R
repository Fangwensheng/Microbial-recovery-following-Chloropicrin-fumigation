# 09_FigureS3.R
source("R/00_setup.R")
a <- read_csv(file.path(D_PANEL,"Supplementary","FigureS3",
                        "FigS3_functional_alpha_raw.csv"),show_col_types=FALSE) %>%
  mutate(Day=Week*7)
pc <- read_csv(file.path(D_PANEL,"Supplementary","FigureS3",
                         "FigS3_functional_PCoA_coordinates.csv"),show_col_types=FALSE)
bs_file <- file.path(D_REC,"cross_domain_bootstrap_summary_Rrerun.csv")
if(!file.exists(bs_file)) bs_file <- file.path(D_REC,"cross_domain_bootstrap_summary_v2.csv")
draws_file <- file.path(D_REC,"cross_domain_bootstrap_draws_Rrerun.csv")
if(!file.exists(draws_file)) draws_file <- file.path(D_REC,"cross_domain_bootstrap_draws_v2.csv")
bs <- read_csv(bs_file,show_col_types=FALSE)
draws <- read_csv(draws_file,show_col_types=FALSE)
# harmonize column names from R rerun vs bundled Python audit result
if("Function_minus_taxonomy" %in% names(bs)) bs <- bs %>% rename(contrast_mean=Function_minus_taxonomy, context=Context, ci_low=CI_low, ci_high=CI_high)

al <- a %>% select(SampleID,Treatment,Day,KO_Shannon,COG_Shannon,KO_Simpson,COG_Simpson) %>%
  pivot_longer(-c(SampleID,Treatment,Day),names_to="Metric",values_to="Value") %>%
  group_by(Treatment,Day,Metric) %>% summarise(mean=mean(Value),se=sd(Value)/sqrt(n()),.groups="drop")
p1 <- ggplot(al,aes(Day,mean,colour=Treatment,linetype=Treatment,shape=Treatment))+
  geom_line(linewidth=.55)+geom_point(size=1.2)+
  facet_wrap(~Metric,scales="free_y",ncol=2)+
  scale_colour_manual(values=treatment_cols)+scale_linetype_manual(values=treatment_ltypes)+
  scale_shape_manual(values=treatment_shapes)+theme_pub()+labs(x="Days after fumigation",y=NULL)

p2 <- ggplot(pc,aes(PCoA1,PCoA2,colour=Treatment,shape=Treatment))+
  geom_path(aes(group=Treatment),alpha=.3)+geom_point(size=1.2)+
  facet_wrap(~Profile,scales="free")+scale_colour_manual(values=treatment_cols)+
  scale_shape_manual(values=treatment_shapes)+theme_pub()+labs(x="PCoA1",y="PCoA2")

# Formal rerun: hierarchical bootstrap slope contrast.
p3 <- ggplot(draws,aes(Slope_contrast_function_minus_taxonomy,fill=Context))+
  geom_density(alpha=.35)+geom_vline(xintercept=0,linetype=2)+
  scale_fill_manual(values=c(Unplanted="#4C6F8A",`Tomato-planted`="#2A9D8F"))+
  theme_pub()+labs(x="Function - taxonomy slope contrast per 10 d",y="Density")
p4 <- ggplot(bs,aes(contrast_mean,context))+
  geom_vline(xintercept=0,linetype=2,colour="grey60")+
  geom_errorbarh(aes(xmin=ci_low,xmax=ci_high),height=.12)+geom_point(size=2)+
  theme_pub()+labs(x="Bootstrap slope contrast (95% CI)",y=NULL)
save_pdf((p1/p2)/(p3|p4)+plot_annotation(tag_levels="a"),
         "FigureS3_R.pdf",width=7.2,height=9.0)
