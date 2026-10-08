# 12_FigureS6.R
source("R/00_setup.R")
long_file <- file.path(D_REC,"provider_replacement_long_Rrerun.csv")
if(file.exists(long_file)) {
  long <- read_csv(long_file,show_col_types=FALSE)
} else {
  long <- read_csv(file.path(D_PANEL,"Supplementary","FigureS6",
                             "FigS6_provider_replacement_long.csv"),show_col_types=FALSE)
}
thr_file <- file.path(D_REC,"provider_replacement_threshold_sensitivity_Rrerun.csv")
if(file.exists(thr_file)){
  sens <- read_csv(thr_file,show_col_types=FALSE) %>%
    transmute(Section="Threshold sensitivity",Context,
              Near_reference_threshold,Percent=Percent_strong,
              N_evaluable,N_strong)
} else {
  sens <- read_csv(file.path(D_PANEL,"Supplementary","FigureS6",
                             "FigS6_threshold_sensitivity_and_summary.csv"),show_col_types=FALSE)
}
rec <- read_csv(file.path(D_RAW,"DESeq2_KO_recovery_descriptors.csv"),show_col_types=FALSE) %>%
  mutate(Context=recode(Mode,Natural="Unplanted",Plant="Tomato-planted"))

# a threshold sensitivity
a0 <- sens %>% filter(Section=="Threshold sensitivity")
a <- ggplot(a0,aes(factor(Near_reference_threshold),Percent,fill=Context))+
  geom_col(position="dodge",width=.65)+
  scale_fill_manual(values=c(Unplanted="#4C6F8A",`Tomato-planted`="#2A9D8F"))+
  theme_pub()+labs(x="Near-reference |log2FC| threshold",
                   y="% near-reference KOs with BC >= 0.50")

# b provider replacement vs reference drift through time
b0 <- long %>% group_by(Context,Day) %>%
  summarise(Treatment_reference=median(Provider_replacement_BC,na.rm=TRUE),
            Reference_drift=median(Reference_drift_from_7d_BC,na.rm=TRUE),.groups="drop") %>%
  pivot_longer(c(Treatment_reference,Reference_drift),
               names_to="Comparison",values_to="BC")
b <- ggplot(b0,aes(Day,BC,colour=Comparison,linetype=Comparison))+
  geom_line(linewidth=.65)+geom_point(size=1.4)+facet_wrap(~Context)+
  theme_pub()+labs(x="Days after fumigation",y="Median Bray-Curtis")

# c integrated provider replacement vs functional residual
trapz <- function(x,y){
  o<-order(x); x<-x[o]; y<-y[o]
  sum(diff(x)*(head(y,-1)+tail(y,-1))/2)/(max(x)-min(x))
}
c0 <- long %>% group_by(KO,Context) %>%
  summarise(Integrated_provider_replacement=trapz(Day,Provider_replacement_BC),.groups="drop") %>%
  left_join(rec %>% select(KO,Context,Time_integrated_absolute_residual_7to70d),
            by=c("KO","Context"))
c <- ggplot(c0,aes(Integrated_provider_replacement,
                   Time_integrated_absolute_residual_7to70d))+
  geom_point(alpha=.18,size=.7)+geom_smooth(method="lm",se=FALSE,colour="black",linewidth=.6)+
  facet_wrap(~Context)+theme_pub()+
  labs(x="Time-integrated provider replacement",y="Time-integrated functional residual")

# d 70-d quadrant sensitivity / reference drift
d <- long %>% filter(Day==70) %>%
  ggplot(aes(Absolute_function_deviation,Provider_replacement_BC,
             colour=Reference_drift_from_7d_BC))+
  geom_vline(xintercept=.25,linetype=2)+geom_hline(yintercept=.5,linetype=2)+
  geom_point(alpha=.5,size=.8)+facet_wrap(~Context)+
  scale_colour_viridis_c()+theme_pub()+
  labs(x=expression("|log"[2]*"FC| at 70 d"),y="Provider replacement",colour="Ref. drift")

save_pdf((a|b)/(c|d)+plot_annotation(tag_levels="a"),"FigureS6_R.pdf",width=7.2,height=6.8)
