# 05_Figure5.R
source("R/00_setup.R")
driver_file <- file.path(D_REC,"Fig5_driver_master_reclassified_Rrerun.csv")
if(!file.exists(driver_file)) driver_file <- file.path(D_REC,"Fig5_driver_master_reclassified_v2.csv")
if(!file.exists(driver_file)) driver_file <- file.path(D_PANEL,"Main_Figures","Figure5","Fig5b_KO_level_driver_raw.csv")
driver <- read_csv(driver_file,show_col_types=FALSE)

univ_file <- file.path(D_REC,"Fig5_univariate_reclassified_Rrerun.csv")
if(!file.exists(univ_file)) univ_file <- file.path(D_REC,"Fig5_univariate_reclassified_v2.csv")
if(!file.exists(univ_file)) univ_file <- file.path(D_PANEL,"Main_Figures","Figure5","Fig5c_univariate_driver_associations.csv")
univ <- read_csv(univ_file,show_col_types=FALSE)
if("Effect" %in% names(univ)) univ <- univ %>% rename(rho=Effect)

partial_file <- file.path(D_REC,"Fig5_partial_spearman_Rrerun.csv")
if(!file.exists(partial_file)) partial_file <- file.path(D_REC,"Fig5_partial_spearman_v2.csv")
if(!file.exists(partial_file)) partial_file <- file.path(D_PANEL,"Main_Figures","Figure5","Fig5c_adjusted_partial_Spearman.csv")
partial <- read_csv(partial_file,show_col_types=FALSE)
if("Effect" %in% names(partial)) partial <- partial %>% rename(partial_rho=Effect)

cv_file <- file.path(D_REC,"Fig5_nestedCV_reclassified_Rrerun.csv")
if(!file.exists(cv_file)) cv_file <- file.path(D_REC,"Fig5_nestedCV_reclassified_v2.csv")
if(!file.exists(cv_file)) cv_file <- file.path(D_PANEL,"Main_Figures","Figure5","Fig5d_nested_CV_ridge.csv")
cv <- read_csv(cv_file,show_col_types=FALSE)
if("Effect" %in% names(cv)) cv <- cv %>% rename(Mean_heldout_R2=Effect)

# a: editable conceptual predictor blocks
boxes <- tibble(
  x=c(1,2,3), label=c("Initial state","Functional architecture","Provider traits"),
  sub=c("day-0 magnitude\ndirection\nbaseline abundance",
        "provider richness\neffective diversity\nphylogenetic breadth\nKO copy",
        "GC / genome size / gene density\nrrn / tRNA / CAZyme\ntransport / stress / repair\nmotility / dormancy")
)
a <- ggplot(boxes,aes(x,1))+
  geom_label(aes(label=paste0(label,"\n",sub)),size=2.5,label.size=.25,
             label.padding=grid::unit(.18,"lines"))+
  annotate("segment",x=1.2,xend=2.8,y=.55,yend=.55,arrow=arrow(length=grid::unit(.12,"cm")))+
  annotate("text",x=2,y=.44,label="Function-specific residual",size=2.7,fontface="bold")+
  coord_cartesian(xlim=c(.5,3.5),ylim=c(.3,1.35),clip="off")+theme_void()

# b1/b2: observed relationships
b1 <- driver %>%
  select(KO,Initial_direction,Unplanted_residual,Tomato_planted_residual) %>%
  pivot_longer(ends_with("_residual"),names_to="Context",values_to="Residual") %>%
  mutate(Context=recode(Context,Unplanted_residual="Unplanted",
                        Tomato_planted_residual="Tomato-planted")) %>%
  ggplot(aes(Initial_direction,Residual,fill=Initial_direction))+
  geom_boxplot(width=.6,outlier.shape=NA)+
  facet_wrap(~Context)+
  scale_fill_manual(values=c(Enriched=COL_FUM,Suppressed=COL_REF))+
  labs(x="Day-0 direction",y="Time-integrated absolute residual")+theme_pub()+
  theme(legend.position="none")

b2 <- driver %>%
  select(Provider_richness,Unplanted_residual,Tomato_planted_residual) %>%
  pivot_longer(ends_with("_residual"),names_to="Context",values_to="Residual") %>%
  mutate(Context=recode(Context,Unplanted_residual="Unplanted",
                        Tomato_planted_residual="Tomato-planted")) %>%
  ggplot(aes(Provider_richness,Residual))+
  geom_point(alpha=.18,size=.8)+
  geom_smooth(method="loess",se=TRUE,linewidth=.6,colour="black")+
  facet_wrap(~Context)+
  labs(x="Provider richness",y="Time-integrated absolute residual")+theme_pub()
b <- b1|b2

# c: forest of univariate Spearman; adjusted redundancy estimates overlaid.
uplot <- univ %>% mutate(Predictor=factor(Predictor,levels=rev(unique(Predictor))))
c <- ggplot(uplot,aes(rho,Predictor,colour=Context,shape=Context))+
  geom_vline(xintercept=0,linetype=2,colour="grey70")+
  geom_point(size=1.7)+
  scale_colour_manual(values=c(Unplanted="#4C6F8A",`Tomato-planted`="#2A9D8F"))+
  labs(x="Spearman \u03c1",y=NULL)+theme_pub()

# d: nested-CV held-out R2
d <- cv %>%
  ggplot(aes(Block,Mean_heldout_R2,fill=Context))+
  geom_hline(yintercept=0,linetype=2,colour="grey60")+
  geom_col(position=position_dodge(.7),width=.62)+
  geom_errorbar(aes(ymin=Mean_heldout_R2-SD_outer_fold_R2,
                    ymax=Mean_heldout_R2+SD_outer_fold_R2),
                position=position_dodge(.7),width=.18,linewidth=.35)+
  scale_fill_manual(values=c(Unplanted="#4C6F8A",`Tomato-planted`="#2A9D8F"))+
  labs(x=NULL,y=expression("Nested-CV held-out "*R^2))+theme_pub()+
  theme(axis.text.x=element_text(angle=25,hjust=1))

p <- a/(b)/(c|d)+plot_annotation(tag_levels="a")
save_pdf(p,"Figure5_R.pdf",width=7.2,height=8.5)
