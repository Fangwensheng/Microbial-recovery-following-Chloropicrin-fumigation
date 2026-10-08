# 04_Figure4.R
source("R/00_setup.R")
h <- read_csv(file.path(D_PANEL,"Main_Figures","Figure4",
                        "Fig4a_oriented_KO_heatmap_data.csv"),show_col_types=FALSE)
traj <- read_csv(file.path(D_PANEL,"Main_Figures","Figure4",
                           "Fig4b_recovery_class_trajectory_summary.csv"),show_col_types=FALSE)
trans <- read_csv(file.path(D_PANEL,"Main_Figures","Figure4",
                            "Fig4c_class_transition_matrix.csv"),show_col_types=FALSE)
res <- read_csv(file.path(D_PANEL,"Main_Figures","Figure4",
                          "Fig4d_time_integrated_residual_pairs.csv"),show_col_types=FALSE)
cls <- c("Rapid","Overshoot/reversal","Gradual/incomplete","Persistent/strengthened")

# If the optional DESeq2 trajectory re-analysis has been run, rebuild Fig.4
# plotting tables directly from the recomputed trajectories/classes/descriptors.
traj_r_file <- file.path(D_REC,"KO_7to70d_oriented_trajectories_Rrerun.csv")
class_r_file <- file.path(D_REC,"KO_trajectory_classes_Rrerun.csv")
desc_r_file <- file.path(D_REC,"KO_recovery_descriptors_Rrerun.csv")
if(file.exists(traj_r_file) && file.exists(class_r_file) && file.exists(desc_r_file)){
  trr <- read_csv(traj_r_file,show_col_types=FALSE)
  clr <- read_csv(class_r_file,show_col_types=FALSE)
  dsr <- read_csv(desc_r_file,show_col_types=FALSE)
  h <- trr %>% mutate(Context=recode(Mode,Natural="Unplanted",Plant="Tomato-planted")) %>%
    left_join(clr %>% select(KO,Natural_class,Plant_class),by="KO") %>%
    mutate(Unplanted_class=Natural_class,Tomato_planted_class=Plant_class) %>%
    select(KO,Context,Day,Oriented_relative_effect,Unplanted_class,Tomato_planted_class)
  traj <- trr %>% mutate(Context=recode(Mode,Natural="Unplanted",Plant="Tomato-planted")) %>%
    left_join(clr %>% select(KO,Natural_class),by="KO") %>%
    group_by(Natural_class,Context,Day) %>%
    summarise(Mean_oriented_effect=mean(Oriented_relative_effect,na.rm=TRUE),
              SE=sd(Oriented_relative_effect,na.rm=TRUE)/sqrt(sum(is.finite(Oriented_relative_effect))),.groups="drop") %>%
    rename(Unplanted_defined_class=Natural_class)
  trans <- clr %>% count(Natural_class,Plant_class,name="N") %>% group_by(Natural_class) %>%
    mutate(Row_percent=100*N/sum(N)) %>% ungroup() %>%
    rename(Unplanted_class=Natural_class,Tomato_planted_class=Plant_class)
  res <- dsr %>% select(KO,Mode,Time_integrated_absolute_residual_7to70d) %>%
    pivot_wider(names_from=Mode,values_from=Time_integrated_absolute_residual_7to70d) %>%
    transmute(KO,Unplanted_integrated_abs_residual=Natural,Tomato_planted_integrated_abs_residual=Plant)
}

# a: two editable heatmaps, KOs ordered by Unplanted class then Unplanted residual-like trajectory.
ord <- h %>% filter(Context=="Unplanted") %>%
  group_by(KO,Unplanted_class) %>% summarise(score=mean(abs(Oriented_relative_effect)),.groups="drop") %>%
  mutate(Unplanted_class=factor(Unplanted_class,levels=cls)) %>%
  arrange(Unplanted_class,score) %>% pull(KO)
h <- h %>% mutate(KO=factor(KO,levels=ord),
                  Context=factor(Context,levels=c("Unplanted","Tomato-planted")))
a <- ggplot(h,aes(Day,KO,fill=Oriented_relative_effect))+
  geom_tile()+facet_grid(.~Context)+
  scale_fill_gradient2(low="#3B6FB6",mid="white",high="#C94C4C",midpoint=0,
                       limits=c(-2,2),oob=squish)+
  scale_x_continuous(breaks=c(7,14,28,42,56,70))+
  labs(x="Days after fumigation",y="1,640 disturbed KOs",fill="Oriented\neffect")+
  theme_pub()+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())

# b class trajectories
b <- traj %>% mutate(Unplanted_defined_class=factor(Unplanted_defined_class,levels=cls)) %>%
  ggplot(aes(Day,Mean_oriented_effect,colour=Context,linetype=Context))+
  geom_hline(yintercept=0,linetype=2,colour="grey60")+
  geom_line(linewidth=.65)+geom_point(size=1.4)+
  geom_errorbar(aes(ymin=Mean_oriented_effect-SE,ymax=Mean_oriented_effect+SE),
                width=2,linewidth=.3)+
  facet_wrap(~Unplanted_defined_class,ncol=2)+
  scale_colour_manual(values=c(Unplanted="#4C6F8A",`Tomato-planted`="#2A9D8F"))+
  labs(x="Days after fumigation",y="Oriented relative effect")+theme_pub()

# c transition matrix with percentages
c <- trans %>% mutate(Unplanted_class=factor(Unplanted_class,levels=cls),
                      Tomato_planted_class=factor(Tomato_planted_class,levels=cls)) %>%
  ggplot(aes(Tomato_planted_class,Unplanted_class,fill=Row_percent))+
  geom_tile(colour="white")+
  geom_text(aes(label=sprintf("%.1f%%",Row_percent)),size=2.2)+
  scale_fill_gradient(low="white",high="#4C78A8")+
  labs(x="Tomato-planted class",y="Unplanted class",fill="Row %")+
  theme_pub()+theme(axis.text.x=element_text(angle=30,hjust=1))

# d paired residual
d <- ggplot(res,aes(Unplanted_integrated_abs_residual,
                    Tomato_planted_integrated_abs_residual))+
  geom_abline(slope=1,intercept=0,linetype=2,colour="grey60")+
  geom_point(alpha=.25,size=1)+
  labs(x="Unplanted cumulative residual",
       y="Tomato-planted cumulative residual")+theme_pub()

p <- a/(b|(c/d))+plot_annotation(tag_levels="a")
save_pdf(p,"Figure4_R.pdf",width=7.2,height=8.8)
